import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
from pysam import VariantFile
from sqlalchemy import func
from sqlmodel import Session, SQLModel, select
from typer.testing import CliRunner

from spliceai_cache.cli import app
from spliceai_cache.db import (
    Contig,
    ContigOrder,
    Position,
    Reference,
    SpliceAIAnnotation,
    SpliceAIConfig,
    SpliceAIInfoField,
    Variant,
)
from spliceai_cache.utilities import SHA256, ContigTuple, sqlite_engine
from spliceai_cache.vcf import extract_vcf, ingest_vcf


@dataclass(frozen=True)
class CacheInputs:
    uri: Path
    dictionary: Path
    annotation: Path
    input: Path
    output: Path

    @property
    def options(self) -> dict[str, Any]:
        return {
            "uri": self.uri,
            "ref_dict": self.dictionary,
            "annotation": self.annotation,
            "distance": 50,
            "mask": False,
        }

    def ingest(self, input_path: Path | None = None, **options: Any) -> tuple[int, int]:
        return ingest_vcf(
            input_path=input_path or self.input,
            **self.options,
            **options,
        )

    def extract(self, **options: Any) -> int:
        return extract_vcf(output_path=self.output, **self.options, **options)


@pytest.fixture
def cache(tmp_path: Path) -> CacheInputs:
    dictionary = tmp_path / "reference.dict"
    dictionary.write_text(
        "@HD\tVN:1.6\n"
        "@SQ\tSN:chr1\tLN:100\tM5:11111111111111111111111111111111\n"
        "@SQ\tSN:chr2\tLN:200\tM5:22222222222222222222222222222222\n"
    )
    annotation = tmp_path / "genes.txt"
    annotation.write_text("# mock SpliceAI annotation source\n")
    vcf = tmp_path / "annotated.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        "##contig=<ID=chr2,length=200>\n"
        '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 scores">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr2\t5\t.\tC\tG\t.\t.\tSpliceAI=G|GENE3|0.1|0|0|0|1|0|0|0\n"
        "chr1\t2\t.\tA\tT,C\t.\t.\tSpliceAI=T|GENE1|0.9|0|0|0|1|0|0|0,C|GENE2|0.2|0|0|0|2|0|0|0,T|GENE4|0.3|0|0|0|3|0|0|0\n"
        "chr1\t3\t.\tG\tA\t.\t.\t.\n"
    )
    inputs = CacheInputs(
        uri=tmp_path / "cache.sqlite",
        dictionary=dictionary,
        annotation=annotation,
        input=vcf,
        output=tmp_path / "extracted.vcf",
    )
    SQLModel.metadata.create_all(sqlite_engine(inputs.uri))
    return inputs


def _count(session: Session, model: type[SQLModel]) -> int:
    return session.exec(select(func.count()).select_from(model)).one()


def test_ingest_is_idempotent_and_extracts_ordered_vcf(cache: CacheInputs) -> None:
    first_config_id, first_count = cache.ingest()
    second_config_id, second_count = cache.ingest()

    assert first_config_id == second_config_id
    assert first_count == second_count == 3
    with Session(sqlite_engine(cache.uri)) as session:
        assert tuple(
            _count(session, model)
            for model in (
                Reference,
                Contig,
                ContigOrder,
                Position,
                Variant,
                SpliceAIInfoField,
                SpliceAIConfig,
                SpliceAIAnnotation,
            )
        ) == (1, 2, 2, 3, 4, 1, 1, 3)

    assert cache.extract() == 3
    with VariantFile(cache.output) as output:
        assert tuple(output.header.contigs) == ("chr1", "chr2")
        records = list(output)
    assert [
        (record.contig, record.pos, record.ref, record.alts) for record in records
    ] == [
        ("chr1", 2, "A", ("C",)),
        ("chr1", 2, "A", ("T",)),
        ("chr2", 5, "C", ("G",)),
    ]
    assert records[0].info["SpliceAI"] == ("C|GENE2|0.2|0|0|0|2|0|0|0",)
    assert records[1].info["SpliceAI"] == (
        "T|GENE1|0.9|0|0|0|1|0|0|0",
        "T|GENE4|0.3|0|0|0|3|0|0|0",
    )


def test_extract_writes_empty_vcf_when_configuration_is_not_cached(
    cache: CacheInputs,
) -> None:
    assert cache.extract() == 0

    with VariantFile(cache.output) as output:
        assert tuple(output.header.contigs) == ("chr1", "chr2")
        assert "SpliceAI" not in output.header.info
        assert list(output) == []


def test_extract_cli_reports_errors_without_a_traceback(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        app,
        [
            "extract",
            "--output",
            str(tmp_path / "output.vcf"),
            "--uri",
            str(tmp_path / "cache.sqlite"),
            "--dict",
            str(tmp_path / "missing.dict"),
            "--annotation",
            str(tmp_path / "missing.txt"),
            "--distance",
            "50",
            "--mask",
        ],
    )

    assert result.exit_code == 1
    assert "Error:" in result.output
    assert "AttributeError" not in result.output


def test_extract_writes_header_without_spliceai_for_empty_config(
    cache: CacheInputs,
) -> None:
    empty_vcf = cache.input.with_name("empty.vcf")
    empty_vcf.write_text(
        "\n".join(
            line for line in cache.input.read_text().splitlines() if line.startswith("#")
        )
        + "\n"
    )
    _config_id, cached_count = cache.ingest(input_path=empty_vcf)

    assert cached_count == 0
    assert cache.extract() == 0
    with VariantFile(cache.output) as output:
        assert tuple(output.header.contigs) == ("chr1", "chr2")
        assert "SpliceAI" not in output.header.info
        assert list(output) == []


def test_all_variants_round_trip_includes_records_without_annotations(
    cache: CacheInputs,
) -> None:
    _config_id, cached_count = cache.ingest(all_variants=True)

    assert cached_count == 4
    with Session(sqlite_engine(cache.uri)) as session:
        assert _count(session, Variant) == 4
        assert _count(session, SpliceAIAnnotation) == 4

    unannotated_vcf = cache.input.with_name("unannotated.vcf")
    lines = cache.input.read_text().splitlines()
    unannotated_vcf.write_text(
        "\n".join(
            line if line.startswith("#") else line.rsplit("\t", 1)[0] + "\t."
            for line in lines
        )
        + "\n"
    )
    _config_id, cached_count = cache.ingest(
        input_path=unannotated_vcf,
        all_variants=True,
    )
    assert cached_count == 4

    assert cache.extract() == 3
    assert cache.extract(all_variants=True) == 4
    with VariantFile(cache.output) as output:
        records = list(output)

    assert [(record.contig, record.pos, record.alts) for record in records] == [
        ("chr1", 2, ("C",)),
        ("chr1", 2, ("T",)),
        ("chr1", 3, ("A",)),
        ("chr2", 5, ("G",)),
    ]
    assert "SpliceAI" not in records[2].info


def test_info_versions_are_stored_separately_and_extracted_together(
    cache: CacheInputs,
) -> None:
    first_id, _cached_count = cache.ingest()
    other_vcf = cache.input.with_name("other-header.vcf")
    other_vcf.write_text(
        cache.input.read_text().replace(
            "SpliceAIv1.3.1 scores", "SpliceAIv1.4.0 scores"
        )
    )
    second_id, _cached_count = cache.ingest(input_path=other_vcf)

    assert first_id != second_id
    with Session(sqlite_engine(cache.uri)) as session:
        assert _count(session, SpliceAIInfoField) == 2
        assert _count(session, SpliceAIConfig) == 2
        versions = session.exec(
            select(SpliceAIInfoField.version).order_by(SpliceAIInfoField.version)
        ).all()
        assert versions == ["1.3.1", "1.4.0"]

    assert cache.extract() == 3
    with VariantFile(cache.output) as output:
        assert output.header.info["SpliceAI"].description == "SpliceAIv1.4.0 scores"


def test_contig_tuple_parses_flexible_sq_tag_order(tmp_path: Path) -> None:
    dictionary = tmp_path / "reference.dict"
    dictionary.write_text(
        "@HD\tVN:1.6\n"
        "@SQ\tLN:10\tM5:AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tSN:chr1\tUR:file.fa\n"
        "@SQ\tSN:chr2\tLN:20\tM5:bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n"
    )

    contigs = ContigTuple.from_path(dictionary)

    assert contigs == (
        {"sn": "chr1", "ln": 10, "md5": "a" * 32},
        {"sn": "chr2", "ln": 20, "md5": "b" * 32},
    )
    assert len(contigs.sha256) == 64
    assert contigs.sha256 == ContigTuple.from_path(dictionary).sha256


def test_contig_tuple_rejects_an_invalid_dictionary(tmp_path: Path) -> None:
    dictionary = tmp_path / "reference.dict"
    dictionary.write_text("@SQ\tSN:chr1\tLN:10\n")

    with pytest.raises(ValueError, match="missing M5"):
        ContigTuple.from_path(dictionary)


def test_sha256_digest_file(tmp_path: Path) -> None:
    path = tmp_path / "annotation.txt"
    path.write_bytes(b"annotation\n")

    assert SHA256.digest_file(path) == hashlib.sha256(b"annotation\n").hexdigest()
