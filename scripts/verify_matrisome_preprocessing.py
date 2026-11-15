"""Trace diagnosis-dependent preprocessing in a pinned public classifier implementation.

Inspect code and the deposited diagnosis/test tables offline. No third-party code is
executed. The receipt describes the public pipeline, not the authors' historical run.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

REVISION = "df538583652d3359d71959437837ea3e0d1195e6"
SOURCES = {
    "data_prep/prep.ipynb": "cc54b17ac7c4c7b1184a84aefa118272461e0fd764b1f66655b76b66432cd5ea",
    "data_prep/batch_correction.r": "20366a374148d732602894f80a53e126224c278dc9492dfbe560d4097efab6be",
    "analysis/regression_classifier.py": "efdb65cb1adc5bc0763420ee71678d8ec8726eb454a6335b61ae55601460349b",
}
TABLES = {
    "coldata.tsv": "a2696216f28f6f80a8e7db716f280c45e4d2e4264730723fb3775b94751ebc4a",
    "test.tsv": "35cbd873dabf6f5567f9a01a2c6fad29850a2920597328edeac762cca1f4f273",
}


def checked(path, pin):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != pin:
        raise ValueError("Source hash mismatch: " + path.name)
    return raw.decode("utf-8")


def lines(text, fragments):
    result = []
    for fragment in fragments:
        hits = [i for i, line in enumerate(text.splitlines(), 1) if fragment in line]
        if len(hits) != 1:
            raise ValueError("Pinned data-flow expression is missing or ambiguous")
        result.append(hits[0])
    if result != sorted(result):
        raise ValueError("Source order changed")
    return result


def audit(source, tables):
    sources = {path: checked(source / path, pin) for path, pin in SOURCES.items()}
    cells = json.loads(sources["data_prep/prep.ipynb"])["cells"]
    steps = {
        "pool_four_studies": (5, ['["GSE4888", "GSE6364", "GSE51981", "GSE29981"]',
                                  "pd.concat(dataframes, axis=1)", "merged_counts_df.to_csv"]),
        "copy_all_labels": (6, ['METADATA_DIR / "coldata.tsv"', "unified_coldata_path.write_text"]),
        "align_expression": (7, ["counts.columns.isin(coldata.index)", "counts.to_csv(unified_counts_path"]),
        "condition_aware_batch_correction": (8, ['DATA_PREP_DIR / "batch_correction.r"',
                                                  "bc_script, str(unified_counts_path), str(unified_coldata_path)",
                                                  "subprocess.run(cmd)"]),
        "create_matrisome_matrix": (9, ["pd.read_csv(unified_counts_path",
                                         "all_matrisome_counts_df.to_csv(all_matrisome_counts_path"]),
        "make_fit_test_tables": (11, ['coldata_df[\'batch\'] == 2', 'coldata_df[\'batch\'] != 2',
                                      'test_dir / f"{phase}_{gene_set}_counts.tsv"',
                                      "sample_counts_df.to_csv(sample_counts_file"]),
    }
    locations = {}
    for step, (index, fragments) in steps.items():
        cell = cells[index]
        if cell["cell_type"] != "code":
            raise ValueError("Expected an executable notebook cell")
        locations[step] = dict(zero_based_cell=index, one_based_lines=lines("".join(cell["source"]), fragments))

    r_lines = lines(sources["data_prep/batch_correction.r"], [
        "expression_file = args[1]", "clinical_file = args[2]", "clinical_df <- read_tsv(clinical_file)",
        "model.matrix(~ condition, data = clinical_df)", "batch <- clinical_df$batch",
        "ComBat(expression_matrix, batch = batch, mod = model_m, ref.batch = 1)",
        "write_tsv(result, file.path(expression_file))"])
    classifier_lines = lines(sources["analysis/regression_classifier.py"], [
        'subdir = "fit" if fit else "test"', 'counts_filename = f"{phase}_{gene_set}_counts.tsv"',
        "for fit in [True, False]", "coldata_path, counts_path = data_paths(phase, gene_set, fit)",
        "genes, joined_df = preprocess_data(counts_path, coldata_path, condition_map)",
        "best_model = optimize_and_train_model", "    evaluate_model(best_model, x_df, y_df"])
    records = {name: list(csv.DictReader(checked(tables / name, pin).splitlines(), delimiter="\t"))
               for name, pin in TABLES.items()}
    all_samples = {r["sample_name"]: r for r in records["coldata.tsv"]}
    held = {k: r for k, r in all_samples.items() if r["batch"] == "2"}
    predictions = {r["sample_name"]: r for r in records["test.tsv"]}
    if len(all_samples) != 216 or len(held) != 37 or held.keys() != predictions.keys():
        raise ValueError("Clinical input and published holdout rows disagree")
    conditions = {"healthy": "0", "endometriosis": "1"}
    common = sum(conditions[r["condition"]] == predictions[k]["true_label"] for k, r in held.items())
    if common != len(held):
        raise ValueError("Condition does not match the evaluation outcome")
    return dict(revision=REVISION, analysis="static source-code and table audit",
                source_sha256=SOURCES, table_sha256=TABLES, notebook_data_flow=locations,
                R_preprocessing_lines=r_lines, python_classifier_lines=classifier_lines,
                clinical_rows_supplied_to_joint_preprocessing=len(all_samples),
                nominal_test_rows=len(held), test_labels_entering_documented_preprocessing=common,
                test_rows_with_diagnosis_held_out_of_documented_preprocessing=0,
                historical_execution_verified=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path, help="Pinned public authors' checkout")
    p.add_argument("tables", type=Path, help="Directory with the audited coldata.tsv and test.tsv")
    args = p.parse_args()
    print(json.dumps(audit(args.source, args.tables), indent=2))
