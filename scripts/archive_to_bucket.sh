#!/usr/bin/env bash
# Encrypt a run directory and upload it to the Hetzner bucket. No secrets in this file:
# expects AGE_RECIPIENT (public key), S3_BUCKET, S3_ENDPOINT and an aws/s3 profile "hetzner" configured from Vault.
set -euo pipefail
run_dir="$1"; stage="${2:-stage0}"; d=$(date +%Y-%m-%d)
name="$(basename "$run_dir")-$d.tar.zst"
tar --zstd -cf - -C "$(dirname "$run_dir")" "$(basename "$run_dir")" | age -r "$AGE_RECIPIENT" > "/tmp/$name.age"
sha256sum "/tmp/$name.age" > "/tmp/$name.age.sha256"
aws --profile hetzner --endpoint-url "${S3_ENDPOINT:?set S3_ENDPOINT}" \
  s3 cp "/tmp/$name.age" "s3://${S3_BUCKET:?set S3_BUCKET}/research/proton-hopping/$stage/$d/"
aws --profile hetzner --endpoint-url "${S3_ENDPOINT:?set S3_ENDPOINT}" \
  s3 cp "/tmp/$name.age.sha256" "s3://${S3_BUCKET:?set S3_BUCKET}/research/proton-hopping/$stage/$d/"
rm -f "/tmp/$name.age" "/tmp/$name.age.sha256"
