# spliceai-cache

`spliceai-cache` stores annotations from SpliceAI output VCFs in SQLite and exports them again without rerunning SpliceAI. Cache configurations are isolated by the reference dictionary, annotation source file, distance, and mask setting.

The reference dictionary must contain `SN`, `LN`, and `M5` tags on every `@SQ` line. The annotation source is the file supplied to SpliceAI with `--annotation`; its SHA-256 digest is used as part of the cache identity.

### Installation
spliceai-cache supports Python 3.10 through 3.13.

The simplest way to install spliceai-cache is through conda:
```sh
conda install -c MOMA-AUH spliceai-cache
```

Alternately, spliceai-cache can be installed from the [github repository](https://github.com/MOMA-AUH/spliceai-cache.git):
```sh
pip install git+https://github.com/MOMA-AUH/spliceai-cache.git
```

## Usage

```console
spliceai-cache init --uri cache.sqlite

spliceai-cache ingest \
  --input spliceai-output.vcf.gz \
  --uri cache.sqlite \
  --dict reference.dict \
  --annotation annotation.txt \
  --distance 50 \
  --no-mask \
  --all-variants

spliceai-cache extract \
  --output cached.vcf.gz \
  --uri cache.sqlite \
  --dict reference.dict \
  --annotation annotation.txt \
  --distance 50 \
  --no-mask \
  --all-variants
```

Use `--mask` instead of `--no-mask` for results produced with SpliceAI masking enabled.

Pass `--all-variants` to both `ingest` and `extract` to retain and export variants that have no `SpliceAI` value. These records are exported without the `SpliceAI` INFO field. Without the option, ingestion and extraction retain their annotation-only behavior.

Re-ingesting the same variants is safe: existing rows are reused and their annotations are updated.
