from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer
from sqlalchemy.exc import SQLAlchemyError
from sqlmodel import SQLModel

from . import db  # noqa: F401
from .utilities import extract_vcf, ingest_vcf, sqlite_engine

app = typer.Typer(
    no_args_is_help=True,
    help="Store and export SpliceAI VCF annotations in a reusable SQLite cache.",
)


@app.command()
def init(
    uri: Annotated[Path, typer.Option("--uri", help="SQLite cache path.")],
) -> None:
    """Create any missing cache tables."""
    SQLModel.metadata.create_all(sqlite_engine(uri))
    typer.echo(f"Initialized cache at {uri}")


@app.command()
def ingest(
    input: Annotated[Path, typer.Option("--input", help="Annotated input VCF.")],
    uri: Annotated[Path, typer.Option("--uri", help="SQLite cache path.")],
    ref_dict: Annotated[
        Path, typer.Option("--dict", help="SAM reference sequence dictionary.")
    ],
    annotation: Annotated[
        Path, typer.Option("--annotation", help="SpliceAI annotation source file.")
    ],
    distance: Annotated[
        int, typer.Option("--distance", min=0, help="SpliceAI maximum distance.")
    ],
    mask: Annotated[
        bool, typer.Option("--mask/--no-mask", help="SpliceAI masking setting.")
    ],
    all_variants: Annotated[
        bool,
        typer.Option(
            "--all-variants",
            help="Cache variants even when they have no SpliceAI annotation.",
        ),
    ] = False,
) -> None:
    """Add the SpliceAI values from an annotated VCF to the cache."""
    try:
        config_id, count = ingest_vcf(
            input_path=input,
            uri=uri,
            ref_dict=ref_dict,
            annotation=annotation,
            distance=distance,
            mask=mask,
            all_variants=all_variants,
        )
    except (OSError, SQLAlchemyError, ValueError, RuntimeError) as exc:
        raise typer.ClickException(str(exc)) from exc
    item_name = "variants" if all_variants else "annotations"
    typer.echo(f"Cached {count} {item_name} (configuration {config_id})")


@app.command()
def extract(
    output: Annotated[Path, typer.Option("--output", help="Output VCF or BCF.")],
    uri: Annotated[Path, typer.Option("--uri", help="SQLite cache path.")],
    ref_dict: Annotated[
        Path, typer.Option("--dict", help="SAM reference sequence dictionary.")
    ],
    annotation: Annotated[
        Path, typer.Option("--annotation", help="SpliceAI annotation source file.")
    ],
    distance: Annotated[
        int, typer.Option("--distance", min=0, help="SpliceAI maximum distance.")
    ],
    mask: Annotated[
        bool, typer.Option("--mask/--no-mask", help="SpliceAI masking setting.")
    ],
    all_variants: Annotated[
        bool,
        typer.Option(
            "--all-variants",
            help="Include cached variants without SpliceAI annotations.",
        ),
    ] = False,
) -> None:
    """Export cached variants"""
    try:
        count = extract_vcf(
            output_path=output,
            uri=uri,
            ref_dict=ref_dict,
            annotation=annotation,
            distance=distance,
            mask=mask,
            all_variants=all_variants,
        )
    except (OSError, SQLAlchemyError, ValueError, RuntimeError) as exc:
        raise typer.ClickException(str(exc)) from exc
    item_name = "variants" if all_variants else "annotations"
    typer.echo(f"Exported {count} {item_name} to {output}")


if __name__ == "__main__":
    app()
