__Data Analysis Steps — Studd data__

These data are on the remote desktop (and Box).

1. Run `compile_odba(filespath,'')` from the remote desktop on the top-most directory (it performs a recursive search for all `.csv` files. This generates `.mat` files with decimated date, ODBA, and temp data.
2. Copy all the generated `.mat` files (e.g., use the filter bar on the top-most directory: `*.mat`); place those files in `Box > Temp` (or separate directory) for analysis.

__Data Analysis Steps — Westrick data__

These data are on the local hard drive (and Box). 2016 files did not contain a Nest classification, so I applied Studd's algorithm in `categorize_odba.m` which writes a new CSV of the same rows, but limited columns to a `nest` folder; these data can then be used with `compile_odba.m`.

1. Run `compile_odba.m` on either the top-level year folder (`ClassifiedData > XXXX`) for 2015 and 2017 or on (`ClassifiedData > XXXX > nest`) for 2016.
2. Then use `copyFiles.m` to transfer those `.mat` files into the temp folder: `copyFiles('...','...','*.mat')` (operates recursively).

__Data Analysis Steps — Dantzer 2019 data__

Rows from the master spreadsheet were associated with the correct files using `fix_2019filenames.m` which exports `2019modexport.csv` with the relevant columns. There is also a commented loop in that file to copy all CSV files to the `2019 Exported Data` folder for easier handling. Since these data came in raw form, they need to be analyzed with nest data using `categorize_odba.m` which will place new CSV files in a `nest` folder. All `*_nest.csv` were compiled once to have access to a MAT-file with compressd rows (in the `T` table); these files are used to get trim dates with 1-minute precision (close enough) without loading the whole original CSV.

1. Files are trimmed using `trim_2019files.m` which saves `Ttrim` in the MAT-file.
2. `Ttrim` is reloaded and used to trim the original CSV files.
3. Trimmed CSV files are categorized into `nest` folder. _Unless new files from 2016 are added, steps 1–3 are done forever._
4. Nest files are compiled to MAT-files that contain T and Tstat.

__Generate Squirrel Key (sqkey)__

The squirrel key places squirrel meta data together and associates them with a filename.

1. Run `build_sqkey.m`
2. Run `view_allSqkey.m`