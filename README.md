__Data Analysis Steps — Studd data__

1. Run `compile_odba(filespath)` from the remote desktop on the top-most directory (it performs a recursive search for all `.csv` files. This generates `.mat` files with decimated date, ODBA, and temp data.
2. Run `run_importKey.m` to add a `filenames` column to the squirrel key. This will save a mat-file named, and with the table, `sqkey` to Box. Move the file to the MATLAB working directory (mat-files are currently ignored).
3. Copy all the generated `.mat` files (e.g., use the filter bar on the top-most directory: `*.mat`); place those files in `Box > Temp` (or separate directory) for analysis.
4. Use `find_dayNight(loadfile,doFig,doWrite)` to generate 'meta' data (e.g., awake, asleep, sunrise, sunset, day length). See `run_dayNightByOdbaPlots.m` for a simple FOR loop.
5. View all data with `analyze_MvsF.m`.