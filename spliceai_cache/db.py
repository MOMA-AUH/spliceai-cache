import base64
import re
from collections.abc import Iterator, Sequence
from pathlib import Path

from sqlalchemy import func
from sqlalchemy.orm import selectinload
from sqlmodel import Field, Relationship, Session, SQLModel, UniqueConstraint, select

from .utilities import SHA256, ContigTuple

SPLICEAI_VERSION_REGEX = re.compile(
    r"\bSpliceAIv(?P<version>\d+(?:\.\d+)*)\b", re.IGNORECASE
)

AnnotationRow = tuple[str, int, str, str, str]


class ContigOrder(SQLModel, table=True):
    __tablename__ = "contig_order"
    __table_args__ = (UniqueConstraint("reference_id", "contig_id"),)

    id: int | None = Field(default=None, primary_key=True, index=True)
    reference_id: int = Field(foreign_key="reference.id", nullable=False, index=True)
    contig_id: int = Field(foreign_key="contig.id", nullable=False, index=True)
    order: int = Field(nullable=False, index=True)

    reference: "Reference" = Relationship(back_populates="contig_orders")
    contig: "Contig" = Relationship(back_populates="contig_orders")


class Contig(SQLModel, table=True):
    __tablename__ = "contig"
    __table_args__ = (UniqueConstraint("sn", "ln", "md5"),)

    id: int | None = Field(default=None, primary_key=True, index=True)
    sn: str = Field(nullable=False, index=True)
    ln: int = Field(nullable=False, index=True)
    md5: str = Field(nullable=False, index=True)

    positions: list["Position"] = Relationship(back_populates="contig")
    contig_orders: list[ContigOrder] = Relationship(back_populates="contig")

    @classmethod
    def populate_from_contig_list(
        cls,
        session: Session,
        contigs: ContigTuple,
    ) -> dict[str, "Contig"]:
        """Return the dictionary contigs, inserting any that are not cached yet."""

        contig_map: dict[str, Contig] = {}
        for contig_data in contigs:
            statement = select(cls).where(
                cls.sn == contig_data["sn"],
                cls.ln == contig_data["ln"],
                cls.md5 == contig_data["md5"],
            )
            contig = session.exec(statement).first()
            if contig is None:
                contig = cls(**contig_data)
                session.add(contig)
                session.flush()
            contig_map[contig.sn] = contig
        return contig_map

    @classmethod
    def populate_from_dict_path(
        cls,
        session: Session,
        dict_path: str | Path,
    ) -> dict[str, "Contig"]:
        """Compatibility wrapper for loading contigs directly from a dictionary."""

        return cls.populate_from_contig_list(session, ContigTuple.from_path(dict_path))

    @classmethod
    def ordered_for_reference(
        cls,
        session: Session,
        reference_id: int,
    ) -> list["Contig"]:
        """Return a reference's contigs in dictionary order."""

        return list(
            session.exec(
                select(cls)
                .join(ContigOrder, ContigOrder.contig_id == cls.id)
                .where(ContigOrder.reference_id == reference_id)
                .order_by(ContigOrder.order)
            ).all()
        )


class Reference(SQLModel, table=True):
    __tablename__ = "reference"

    id: int | None = Field(default=None, primary_key=True, index=True)
    sha256: str = Field(unique=True, nullable=False, index=True)

    contig_orders: list[ContigOrder] = Relationship(
        back_populates="reference",
        sa_relationship_kwargs={
            "cascade": "all, delete-orphan",
            "order_by": "ContigOrder.order",
        },
    )
    spliceai_configs: list["SpliceAIConfig"] = Relationship(back_populates="reference")

    @classmethod
    def find(
        cls,
        session: Session,
        source: str | Path | ContigTuple,
    ) -> "Reference | None":
        """Find the cached reference represented by a sequence dictionary."""

        contigs = (
            source if isinstance(source, ContigTuple) else ContigTuple.from_path(source)
        )
        return session.exec(select(cls).where(cls.sha256 == contigs.sha256)).first()

    @classmethod
    def resolve(
        cls,
        session: Session,
        source: str | Path | ContigTuple,
    ) -> tuple["Reference", dict[str, Contig]]:
        """Resolve a reference and its name-to-contig mapping in one transaction."""

        contigs = (
            source if isinstance(source, ContigTuple) else ContigTuple.from_path(source)
        )
        contig_map = Contig.populate_from_contig_list(session, contigs)
        reference = cls.find(session, contigs)
        if reference is None:
            reference = cls(
                sha256=contigs.sha256,
                contig_orders=[
                    ContigOrder(contig=contig_map[contig["sn"]], order=order)
                    for order, contig in enumerate(contigs, start=1)
                ],
            )
            session.add(reference)
            session.flush()
        else:
            existing = [
                (contig_order.contig_id, contig_order.order)
                for contig_order in reference.contig_orders
            ]
            expected = [
                (contig_map[contig["sn"]].id, order)
                for order, contig in enumerate(contigs, start=1)
            ]
            if existing != expected:
                raise ValueError(
                    "Reference hash matched but its stored contig order is inconsistent"
                )
        return reference, contig_map

    @classmethod
    def get_or_create(cls, session: Session, dict_path: str | Path) -> "Reference":
        """Return the reference represented by a SAM sequence dictionary."""

        reference, _contig_map = cls.resolve(session, dict_path)
        return reference


class Position(SQLModel, table=True):
    __tablename__ = "position"
    __table_args__ = (UniqueConstraint("contig_id", "pos"),)

    id: int | None = Field(default=None, primary_key=True, index=True)
    contig_id: int = Field(foreign_key="contig.id", nullable=False, index=True)
    pos: int = Field(nullable=False)
    ref: str = Field(nullable=False)

    contig: Contig = Relationship(back_populates="positions")
    variants: list["Variant"] = Relationship(back_populates="position")

    @classmethod
    def get_or_create(
        cls,
        session: Session,
        contig: Contig,
        pos: int,
        ref: str,
    ) -> "Position":
        statement = select(cls).where(cls.contig_id == contig.id, cls.pos == pos)
        position = session.exec(statement).first()
        if position is None:
            position = cls(contig=contig, pos=pos, ref=ref)
            session.add(position)
            session.flush()
        elif position.ref != ref:
            raise ValueError(
                f"Conflicting reference alleles at {contig.sn}:{pos}: "
                f"cached {position.ref!r}, received {ref!r}"
            )
        return position


class Variant(SQLModel, table=True):
    __tablename__ = "variant"
    __table_args__ = (UniqueConstraint("pos_id", "alt"),)

    id: int | None = Field(default=None, primary_key=True, index=True)
    pos_id: int = Field(foreign_key="position.id", nullable=False, index=True)
    alt: str = Field(nullable=False)

    position: Position = Relationship(back_populates="variants")
    spliceai_annotations: list["SpliceAIAnnotation"] = Relationship(
        back_populates="variant"
    )

    @classmethod
    def get_or_create(
        cls,
        session: Session,
        position: Position,
        alt: str,
    ) -> "Variant":
        statement = select(cls).where(cls.pos_id == position.id, cls.alt == alt)
        variant = session.exec(statement).first()
        if variant is None:
            variant = cls(position=position, alt=alt)
            session.add(variant)
            session.flush()
        return variant


class SpliceAIInfoField(SQLModel, table=True):
    __tablename__ = "spliceai_info_field"

    id: int | None = Field(default=None, primary_key=True, index=True)
    version: str = Field(nullable=False)
    info_field_base64: str = Field(nullable=False)

    configs: list["SpliceAIConfig"] = Relationship(back_populates="info_field")

    @property
    def header_line(self) -> str:
        """Decode and validate the stored VCF INFO definition."""

        try:
            line = base64.b64decode(self.info_field_base64, validate=True).decode(
                "utf-8"
            )
        except (ValueError, UnicodeDecodeError) as exc:
            raise ValueError(
                "Cached SpliceAI INFO definition is not valid base64"
            ) from exc
        if not line.startswith("##INFO=<"):
            raise ValueError("Cached SpliceAI INFO definition is invalid")
        return line

    @property
    def version_key(self) -> tuple[int, ...]:
        """Return a key that sorts numeric versions after unknown versions."""

        if re.fullmatch(r"\d+(?:\.\d+)*", self.version) is None:
            return ()
        return tuple(int(part) for part in self.version.split("."))

    @classmethod
    def from_header_line(cls, session: Session, line: str) -> "SpliceAIInfoField":
        """Return the cached representation of a VCF SpliceAI INFO definition."""

        match = SPLICEAI_VERSION_REGEX.search(line)
        version = match.group("version") if match else "unknown"
        encoded_line = base64.b64encode(line.encode("utf-8")).decode("ascii")
        return cls.get_or_create(
            session,
            version=version,
            info_field_base64=encoded_line,
        )

    @classmethod
    def get_or_create(
        cls,
        session: Session,
        version: str,
        info_field_base64: str,
    ) -> "SpliceAIInfoField":
        instance = session.exec(
            select(cls).where(
                cls.version == version,
                cls.info_field_base64 == info_field_base64,
            )
        ).first()
        if instance is None:
            instance = cls(version=version, info_field_base64=info_field_base64)
            session.add(instance)
            session.flush()
        return instance


class SpliceAIConfig(SQLModel, table=True):
    __tablename__ = "spliceai_config"
    __table_args__ = (
        UniqueConstraint(
            "reference_id",
            "info_field_id",
            "annotation_sha256",
            "distance",
            "mask",
        ),
    )

    id: int | None = Field(default=None, primary_key=True, index=True)
    reference_id: int = Field(foreign_key="reference.id", nullable=False, index=True)
    info_field_id: int = Field(
        foreign_key="spliceai_info_field.id", nullable=False, index=True
    )
    annotation_sha256: str = Field(nullable=False)
    distance: int = Field(nullable=False)
    mask: bool = Field(nullable=False)

    reference: Reference = Relationship(back_populates="spliceai_configs")
    info_field: SpliceAIInfoField = Relationship(back_populates="configs")
    spliceai_annotations: list["SpliceAIAnnotation"] = Relationship(
        back_populates="config"
    )

    @classmethod
    def matching_inputs(
        cls,
        session: Session,
        ref_dict: str | Path,
        annotation: str | Path,
        distance: int,
        mask: bool,
    ) -> list["SpliceAIConfig"]:
        """Return configurations matching the supplied SpliceAI inputs."""

        reference = Reference.find(session, ref_dict)
        if reference is None or reference.id is None:
            raise ValueError("Reference dictionary is not present in the cache")

        configs = session.exec(
            select(cls)
            .options(selectinload(cls.info_field))
            .where(
                cls.reference_id == reference.id,
                cls.annotation_sha256 == SHA256.digest_file(annotation),
                cls.distance == distance,
                cls.mask == mask,
            )
            .order_by(cls.id)
        ).all()
        if not configs:
            raise ValueError("No cache configuration matches the supplied inputs")
        return list(configs)

    @classmethod
    def newest_info_config(
        cls, configs: Sequence["SpliceAIConfig"]
    ) -> "SpliceAIConfig":
        """Return the configuration with the newest SpliceAI INFO definition."""

        if not configs:
            raise ValueError("At least one cache configuration is required")
        return max(
            configs,
            key=lambda config: config.info_field.version_key,
        )

    @classmethod
    def get_or_create(
        cls,
        session: Session,
        reference_id: int,
        info_field_id: int,
        annotation_sha256: str,
        distance: int,
        mask: bool,
    ) -> "SpliceAIConfig":
        statement = (
            select(cls)
            .where(
                cls.reference_id == reference_id,
                cls.info_field_id == info_field_id,
                cls.annotation_sha256 == annotation_sha256,
                cls.distance == distance,
                cls.mask == mask,
            )
            .order_by(cls.id)
        )
        instance = session.exec(statement).first()
        if instance is None:
            instance = cls(
                reference_id=reference_id,
                info_field_id=info_field_id,
                annotation_sha256=annotation_sha256,
                distance=distance,
                mask=mask,
            )
            session.add(instance)
            session.flush()
        return instance

    @classmethod
    def get_or_create_for_inputs(
        cls,
        session: Session,
        reference: Reference,
        info_field: SpliceAIInfoField,
        annotation: str | Path,
        distance: int,
        mask: bool,
    ) -> "SpliceAIConfig":
        """Resolve a configuration from its input models and annotation file."""

        if reference.id is None:
            raise RuntimeError("Reference did not receive a database id")
        if info_field.id is None:
            raise RuntimeError("SpliceAI INFO field did not receive a database id")
        return cls.get_or_create(
            session,
            reference_id=reference.id,
            info_field_id=info_field.id,
            annotation_sha256=SHA256.digest_file(annotation),
            distance=distance,
            mask=mask,
        )


class SpliceAIAnnotation(SQLModel, table=True):
    """A configuration-specific result; an empty value marks an unannotated variant."""

    __tablename__ = "spliceai_annotation"
    __table_args__ = (UniqueConstraint("variant_id", "config_id"),)

    id: int | None = Field(default=None, primary_key=True, index=True)
    variant_id: int = Field(foreign_key="variant.id", nullable=False, index=True)
    config_id: int = Field(foreign_key="spliceai_config.id", nullable=False, index=True)
    spliceai: str = Field(nullable=False)

    variant: Variant = Relationship(back_populates="spliceai_annotations")
    config: SpliceAIConfig = Relationship(back_populates="spliceai_annotations")

    @classmethod
    def iter_for_configs(
        cls,
        session: Session,
        configs: Sequence[SpliceAIConfig],
        include_unannotated: bool = False,
    ) -> Iterator[AnnotationRow]:
        """Yield distinct variant annotations in reference dictionary order."""

        if not configs:
            raise ValueError("At least one cache configuration is required")
        config_ids = [config.id for config in configs if config.id is not None]
        if len(config_ids) != len(configs):
            raise RuntimeError("Cached configuration has no database id")
        reference_ids = {config.reference_id for config in configs}
        if len(reference_ids) != 1:
            raise ValueError("Cache configurations belong to different references")

        statement = (
            select(
                Contig.sn,
                Position.pos,
                Position.ref,
                Variant.alt,
                func.max(cls.spliceai).label("spliceai"),
            )
            .join(ContigOrder, ContigOrder.contig_id == Contig.id)
            .join(Position, Position.contig_id == Contig.id)
            .join(Variant, Variant.pos_id == Position.id)
            .join(cls, cls.variant_id == Variant.id)
            .where(
                ContigOrder.reference_id == reference_ids.pop(),
                cls.config_id.in_(config_ids),
            )
        )
        if not include_unannotated:
            statement = statement.where(cls.spliceai != "")
        statement = statement.group_by(
            ContigOrder.order,
            Contig.sn,
            Position.pos,
            Position.ref,
            Variant.alt,
        ).order_by(ContigOrder.order, Position.pos, Variant.alt)

        yield from session.exec(statement).yield_per(1000)

    @classmethod
    def upsert(
        cls,
        session: Session,
        variant: Variant,
        config: SpliceAIConfig,
        spliceai: str,
    ) -> "SpliceAIAnnotation":
        statement = select(cls).where(
            cls.variant_id == variant.id, cls.config_id == config.id
        )
        annotation = session.exec(statement).first()
        if annotation is None:
            annotation = cls(variant=variant, config=config, spliceai=spliceai)
            session.add(annotation)
        elif spliceai or not annotation.spliceai:
            annotation.spliceai = spliceai
            session.add(annotation)
        return annotation
