"""Compare saved TSV runs across sessions, without launching a solver."""
import argparse
import json
from pathlib import Path

from run_tsv_validation import ROOT, compare, save, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=ROOT / "gpu/validation/tsv_copper_t40_20261029")
    parser.add_argument("--current", type=Path, default=ROOT / "gpu/validation/tsv_copper_idle_t40_20261029")
    args = parser.parse_args()
    old = json.loads((args.reference / "summary.json").read_text(encoding="utf-8"))
    new = json.loads((args.current / "summary.json").read_text(encoding="utf-8"))
    identities = [json.loads((p / "source_identity.json").read_text(encoding="utf-8")) for p in [args.reference, args.current]]
    assert identities[0] == identities[1], "Source/binary identity changed across runs"
    by_case = {r["case"]: r for r in old["cases"]}
    results = []
    for current in new["cases"]:
        case = current["case"]
        reference = by_case[case]
        record = dict(case=case, backends={})
        for backend in ["cpu", "cuda"]:
            manifests = [json.loads((p / case / backend / "input_manifest.json").read_text(encoding="utf-8"))
                         for p in [args.reference, args.current]]
            for key in ["copied_card_sha256", "input_sha256", "threads", "global_material", "steps", "tolerance"]:
                assert manifests[0][key] == manifests[1][key], (case, backend, key)
            a, b = reference["runs"][backend], current["runs"][backend]
            assert a["success"] and b["success"], "Only compare completed successful runs"
            hashes = [sha256(Path(run["state_file"])) for run in [a, b]]
            if hashes[0] == hashes[1]:
                comparison = dict(rows=current["elements"] * 1001, exactly_equal=True,
                                  comparison_method="Identical full state CSV SHA256", state_sha256=hashes[0],
                                  fields={field: dict(relative_difference=0., max_absolute_difference=0.) for field in ["u", "t"]})
            else:
                comparison, per_step = compare(a, b, current["elements"])
                comparison["exactly_equal"] = all(v["max_absolute_difference"] == 0 for v in comparison["fields"].values())
                per_step.to_csv(args.current / case / ("repeat_" + backend + "_comparison.csv"), index=False)
            record["backends"][backend] = comparison
            print(f"{case} {backend}: exact={comparison['exactly_equal']} u={comparison['fields']['u']['relative_difference']:.6e} t={comparison['fields']['t']['relative_difference']:.6e}", flush=True)
        results.append(record)
        save(args.current / "repeat_comparison.json", dict(reference=str(args.reference), current=str(args.current),
             source_binary_identity_equal=True, copied_input_identity_equal=True, cases=results))


if __name__ == "__main__":
    main()
