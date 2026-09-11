from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from pysam import VariantFile, VariantHeader
from sqlmodel import Session

from .utilities import ContigTuple, sqlite_engine

SPLICEAI_INFO_FIELD = "SpliceAI"


def _info_header_line(vcf: VariantFile) -> str:
    for record in vcf.header.records:
        if record.type == "INFO" and record.get("ID") == SPLICEAI_INFO_FIELD:
            return str(record).rstrip("\r\n")
    raise ValueError(f"Input VCF has no {SPLICEAI_INFO_FIELD!r} INFO definition")


def _spliceai_values(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,)
    if isinstance(value, Iterable):
        return tuple(str(item) for item in value)
    return (str(value),)


def ingest_vcf(
    input_path: str | Path,
    uri: str | Path,
    ref_dict: str | Path,
    annotation: str | Path,
    distance: int,
    mask: bool,
    all_variants: bool = False,
) -> tuple[int, int]:
    from .db import (
        Position,
        Reference,
        SpliceAIAnnotation,
        SpliceAIConfig,
        SpliceAIInfoField,
        Variant,
    )

    if distance < 0:
        raise ValueError("distance must be non-negative")

    contigs = ContigTuple.from_path(ref_dict)
    engine = sqlite_engine(uri)

    cached_count = 0
    with Session(engine) as session, VariantFile(str(input_path)) as vcf:
        info_header_line = _info_header_line(vcf)
        info_field = SpliceAIInfoField.from_header_line(session, info_header_line)
        reference, contig_map = Reference.resolve(session, contigs)
        config = SpliceAIConfig.get_or_create_for_inputs(
            session,
            reference=reference,
            info_field=info_field,
            annotation=annotation,
            distance=distance,
            mask=mask,
        )

        for record in vcf:
            contig = contig_map.get(record.contig)
            if contig is None:
                raise ValueError(
                    f"VCF contig {record.contig!r} is absent from the reference dictionary"
                )
            if record.ref is None:
                raise ValueError(
                    f"VCF record at {record.contig}:{record.pos} has no REF"
                )

            values = _spliceai_values(record.info.get(SPLICEAI_INFO_FIELD))
            values_by_alt: dict[str, list[str]] = {}
            for value in values:
                annotated_alt, separator, _rest = value.partition("|")
                if separator:
                    values_by_alt.setdefault(annotated_alt, []).append(value)

            position = Position.get_or_create(
                session, contig=contig, pos=record.pos, ref=record.ref
            )
            for alt in record.alts or ():
                variant = Variant.get_or_create(session, position=position, alt=alt)
                matching_values = values_by_alt.get(alt)
                if not matching_values and not all_variants:
                    continue
                SpliceAIAnnotation.upsert(
                    session,
                    variant=variant,
                    config=config,
                    spliceai=",".join(matching_values or ()),
                )
                cached_count += 1

        session.commit()
        if config.id is None:
            raise RuntimeError("SpliceAI configuration did not receive a database id")
        return config.id, cached_count


def _output_mode(output_path: str | Path) -> str:
    name = str(output_path).lower()
    if name.endswith(".bcf"):
        return "wb"
    if name.endswith((".vcf.gz", ".vcf.bgz")):
        return "wz"
    return "w"


def _output_header(session: Session, config: object) -> VariantHeader:
    from .db import Contig

    header = VariantHeader()
    for contig in Contig.ordered_for_reference(
        session, reference_id=config.reference_id
    ):
        header.contigs.add(contig.sn, length=contig.ln, md5=contig.md5)
    header.add_line(config.info_field.header_line)
    header.add_meta("source", value="spliceai-cache")
    return header


def extract_vcf(
    output_path: str | Path,
    uri: str | Path,
    ref_dict: str | Path,
    annotation: str | Path,
    distance: int,
    mask: bool,
    all_variants: bool = False,
) -> int:
    from .db import SpliceAIAnnotation, SpliceAIConfig

    if distance < 0:
        raise ValueError("distance must be non-negative")

    engine = sqlite_engine(uri)
    with Session(engine) as session:
        configs = SpliceAIConfig.matching_inputs(
            session,
            ref_dict=ref_dict,
            annotation=annotation,
            distance=distance,
            mask=mask,
        )
        header_config = SpliceAIConfig.newest_info_config(configs)
        header = _output_header(session, header_config)
        rows = SpliceAIAnnotation.iter_for_configs(
            session,
            configs=configs,
            include_unannotated=all_variants,
        )

        output_count = 0
        with VariantFile(
            str(output_path), mode=_output_mode(output_path), header=header
        ) as output_vcf:
            for contig, pos, ref, alt, spliceai in rows:
                output_record = output_vcf.new_record(
                    contig=contig,
                    start=pos - 1,
                    alleles=(ref, alt),
                )
                if spliceai:
                    output_record.info[SPLICEAI_INFO_FIELD] = tuple(spliceai.split(","))
                output_vcf.write(output_record)
                output_count += 1
        return output_count
