from __future__ import annotations

import hashlib
import re
from functools import cached_property
from pathlib import Path
from typing import TypedDict

from sqlalchemy import event
from sqlalchemy.engine import Engine
from sqlmodel import create_engine

SPLICEAI_INFO_FIELD = "SpliceAI"
MD5_REGEX = re.compile(r"^[a-fA-F0-9]{32}$")


class ContigDict(TypedDict):
    """The reference metadata retained for one SAM sequence record."""

    sn: str
    ln: int
    md5: str


class ContigTuple(tuple[ContigDict, ...]):
    """An ordered, canonical representation of the contigs in a SAM dictionary."""

    @classmethod
    def from_path(cls, dict_path: str | Path) -> ContigTuple:
        contigs: list[ContigDict] = []
        seen_names: set[str] = set()

        with Path(dict_path).open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\r\n").split("\t")
                if not fields or fields[0] != "@SQ":
                    continue

                tags: dict[str, str] = {}
                for field in fields[1:]:
                    tag, separator, value = field.partition(":")
                    if separator:
                        tags[tag] = value

                missing = {"SN", "LN", "M5"} - tags.keys()
                if missing:
                    missing_text = ", ".join(sorted(missing))
                    raise ValueError(
                        f"Invalid @SQ record on line {line_number}: missing {missing_text}"
                    )

                name = tags["SN"]
                if not name:
                    raise ValueError(
                        f"Invalid @SQ record on line {line_number}: SN is empty"
                    )
                if name in seen_names:
                    raise ValueError(f"Duplicate contig name in dictionary: {name}")

                try:
                    length = int(tags["LN"])
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid @SQ record on line {line_number}: LN is not an integer"
                    ) from exc
                if length <= 0:
                    raise ValueError(
                        f"Invalid @SQ record on line {line_number}: LN must be positive"
                    )

                md5 = tags["M5"].lower()
                if MD5_REGEX.fullmatch(md5) is None:
                    raise ValueError(
                        f"Invalid @SQ record on line {line_number}: M5 is not an MD5 digest"
                    )

                contigs.append({"sn": name, "ln": length, "md5": md5})
                seen_names.add(name)

        if not contigs:
            raise ValueError(
                f"No @SQ records found in reference dictionary: {dict_path}"
            )
        return cls(contigs)

    @cached_property
    def sha256(self) -> str:
        """Hash the ordered semantic content rather than incidental file formatting."""

        digest = hashlib.sha256()
        for contig in self:
            for value in (contig["sn"], str(contig["ln"]), contig["md5"]):
                encoded = value.encode("utf-8")
                digest.update(len(encoded).to_bytes(8, byteorder="big"))
                digest.update(encoded)
        return digest.hexdigest()


class SHA256:
    @staticmethod
    def digest_file(path: str | Path, block_size: int = 65536) -> str:
        sha = hashlib.sha256()
        with Path(path).open("rb") as handle:
            while block := handle.read(block_size):
                sha.update(block)
        return sha.hexdigest()


def sqlite_engine(uri: str | Path) -> Engine:
    """Create the cache engine and enforce declared foreign keys."""

    engine = create_engine(f"sqlite:///{Path(uri)}")

    @event.listens_for(engine, "connect")
    def enable_sqlite_foreign_keys(dbapi_connection, _connection_record) -> None:
        cursor = dbapi_connection.cursor()
        cursor.execute("PRAGMA foreign_keys=ON")
        cursor.close()

    return engine
