# Experiment and production scripts (as run)

These are the scripts that produced the published results, kept as run so that the numbers in
`docs/VALIDATION.md` can be traced. They expect the Colab layout (Drive mounted, data under `$HRVPP_ROOT`,
default `/content/drive/MyDrive/HRVPP_T34UGV_recovery`) and, where Earth Engine is used, `$EE_PROJECT`.
For new runs use the command-line scripts in [`scripts/`](../../scripts), which wrap the same library code.

| script | what it did | results |
|---|---|---|
| `soil_dvi_sensitivity.py` | iteration 2: soil DVI 0.05–0.11 on 4 × 4 km blocks in the overlap strips | `reports/08_*` |
| `production_learned_correction.py` | iteration 3: learned correction per year on the four T34UGV overlap strips, T34UGV 2017–2024 and the extension 2017 | `reports/11_*` |
| `production_local_gpu.py` | reconstruction of the extension without Earth Engine (Planetary Computer, Colab GPU), 2021 and 2024 | `reports/12_*` |
| `production_colleague_ee.py` | the extension 2018–2024 from Earth Engine exports run in the Code Editor under another account's quota | `reports/13_*` |
