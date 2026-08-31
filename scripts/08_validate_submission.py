from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

EXPECTED_COLUMNS = [
    "proband_id",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "epcr",
    "finding_type",
    "notes",
]

ALLELE_PATTERN = re.compile(r"^[ACGTN]+$")
CHROM_PATTERN = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$")


def validate(path: Path) -> list[str]:
    errors: list[str] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != EXPECTED_COLUMNS:
            errors.append(
                f"Header mismatch.\nExpected: {EXPECTED_COLUMNS}\nActual:   {reader.fieldnames}"
            )
            return errors
        rows = list(reader)

    if not rows:
        errors.append("Submission has no prediction rows.")
    if len(rows) > 10:
        errors.append(f"Submission has {len(rows)} rows; maximum is 10.")

    previous_epcr = float("inf")
    seen_hypotheses: set[tuple[str, ...]] = set()

    for number, row in enumerate(rows, start=2):
        prefix = f"CSV line {number}"

        if row["proband_id"] != "PROBAND01":
            errors.append(f"{prefix}: proband_id must be PROBAND01.")

        try:
            epcr = float(row["epcr"])
            if not 0 < epcr <= 1:
                errors.append(f"{prefix}: epcr must be in (0, 1].")
            if epcr >= previous_epcr:
                errors.append(f"{prefix}: epcr values must be strictly descending.")
            previous_epcr = epcr
        except ValueError:
            errors.append(f"{prefix}: invalid epcr '{row['epcr']}'.")

        if row["finding_type"] not in {"primary", "secondary"}:
            errors.append(f"{prefix}: finding_type must be primary or secondary.")

        for suffix in ("1", "2"):
            chrom = row[f"chrom_{suffix}"]
            pos = row[f"pos_{suffix}"]
            ref = row[f"ref_{suffix}"]
            alt = row[f"alt_{suffix}"]
            values = (chrom, pos, ref, alt)

            if suffix == "1" and not all(values):
                errors.append(f"{prefix}: first variant is incomplete.")
                continue
            if suffix == "2" and any(values) and not all(values):
                errors.append(f"{prefix}: second variant must be fully populated or blank.")
                continue
            if not any(values):
                continue

            if not CHROM_PATTERN.match(chrom):
                errors.append(f"{prefix}: invalid chromosome '{chrom}'.")
            try:
                if int(pos) <= 0:
                    raise ValueError
            except ValueError:
                errors.append(f"{prefix}: invalid position '{pos}'.")
            if not ALLELE_PATTERN.match(ref):
                errors.append(f"{prefix}: invalid REF allele '{ref}'.")
            if not ALLELE_PATTERN.match(alt):
                errors.append(f"{prefix}: invalid ALT allele '{alt}'.")
            if ref == alt:
                errors.append(f"{prefix}: REF and ALT are identical.")

        hypothesis = tuple(row[column] for column in EXPECTED_COLUMNS[1:9])
        if hypothesis in seen_hypotheses:
            errors.append(f"{prefix}: duplicate variant hypothesis.")
        seen_hypotheses.add(hypothesis)

    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("submission", type=Path)
    args = parser.parse_args()

    errors = validate(args.submission)
    if errors:
        print("Submission validation failed:")
        for error in errors:
            print(f"- {error}")
        raise SystemExit(1)

    with args.submission.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    print(f"Submission is valid: {len(rows)} rows, descending EPCR, official schema.")


if __name__ == "__main__":
    main()
