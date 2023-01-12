# MSEED Plot Tool
  - [1 - System requirements:](#1---system-requirements)
  - [2 - Basic usage](#2---basic-usage)
  - [3 - Python functions](#3---python-functions)
    - [3.1 - Main](#31---main)
    - [3.2 - Info](#32---info)
    - [Files management](#files-management)
    - [3.3 - Metadata gathering](#33---metadata-gathering)
    - [3.4 - Elaboration and export](#34---elaboration-and-export)
<br>

## 1 - System requirements
- python3
- numpy
- scipy
- matplotlib >= 3.1.1
- obspy >= 1.1.1


## 2 - Basic usage
The script is used to analyze and plot one or more **miniSEED** files. To get the script to work properly you will need to set some options explained right below:

- `--help`: This option tells the script to print an help message and exit

- `--version`: Prints script's version.

- `--print_output`: Enables script's output on video.

Following options take 1 argument to work correctly:

- `--input`: This option is used to specify a path where the script will search for files to read. You can also use **wildcards** to refine your search. If files like `LOG`, `RNG` and `SOH` are present, they will be discarded as these do not contain trace data to be analyzed.

- `--output`: This option is used to specify a path where the script will save exported data. If the specified directory does not exist, the script will try to create it alog with its subpaths where different kind of export will be stored

    **NOTE:** if this option is not specified, the script will create a directory named "Exports" in the same directory where the script has been placed.

- `--stationxml`: This should be an **XML** created using the FDSN StationXML Schema. This file  contains the station's topology and the instruments response informations. It will be used to gather informations that are necessary to elaborate traces.

- `--archive`  This option enables *archive mode*, it uses a *tar.gz* archive (to be passed as argument) which contains traces data and informations about the datalogger used to export those data. If this option is used the script will automatically calculate the datalogger's frequency response and will elaborate data according to the export and plot options specified. **NOTE:** it this option is use the `--input` e `--stationxml` options will be ignored.

- `--export`: This one will set the output for decovolution operation and PPSD calculation. It can be one of the following
  - `ACC` - Acceleration: output will be in **m/s**<sup>**2**</sup>
  - `VEL` - Velocity: output will be in **m/s**
  - `DISP` - Displacement: output will be in **m**
  
  **NOTE:** if this option is not set it will default to `ACC`

- `--start_time` and `--end_time`: These options can be set to specify a **time window** in which all the traces will be sliced. If not set all the traces will sliced to start and end all the same time: *latest start* and *earliest end*

- `--iris`: This option enables the user to specify an alternative path for the **NRL IRIS Catalog** , by default the script looks for a directory named **"IRIS"** in the same directory where the script has been placed, but you can change this behaviour using this option.

    **NOTE:** Time string should have the following format: `2019-10-29T15:10:30`

Following options are used ad flags, their presence enables the associated kind of plot or modify the script behaviour

- `--plot`, `--plot_raw`, `--plot_ppsd` and `--plot_psd`: These options will tell the script which kind of plot should be exported. They work as follows:
  - `--plot`: This will plot **deconvoluted traces** obtained by removing the instrument's response through the StationXML file.
  - `--plot_raw`: This will plot **raw traces** as they are in the mseed file, with no elaboration.
  - `--plot_psd`: This will plot **PSD**s (Power Spectral Density) of the traces using *Welch's* method provided by the *scipy* library.
  - `--ascii`: This will export deconvoluted traces into an *ASCII* file with parallel columns.
  - `--plot_ppsd`: This will plot **PPSD**s (Probabilistic Power Spectral Density) of the traces using *StationXML* file's informations.

-  `--structural`: If specified, this option enables **PPSD** (Probabilistic Power Spectral Density) plot and sets it to be run with parameters that are suitable for **structural analysis**, otherwise, if the plot will be performed, the parameters will be set for **seismic analysis**.

- `--ascii`: This option enables the export of deconvoluted traces in ASCII format.

- `--psd_ascii`: This option enables the export of **PSD** (Power Spectral Density) in ASCII format.

- `--dayplot`: This option tells the script to plot raw traces in daily format.

- `--spectrogram`: This option tells the script to plot spectrograms for the traces.

- `--seismic`: If specified, this option enables **PSD** plot and sets it to be run with parameters that are suitable for **seismic analysis**, otherwise, if the plot will be performed, the parameters will be set for **structural analysis**.

- `--sds`: If specified, this option tells the script to treat the **input path** as an SDS directory and will search traces data according to that standard

- `--ascii-semicolon`: If present, the semicolon won't be used as end-of-line character in ASCII eport files

Following options take more than 1 argument to work correctly:

- `--unit <label> <conv_factor>`: If present, traces will be converted by multiplying values for `<factor>` value and some plots/exports will use converted traces applying `<label>` as unit of measure.

- `--prefilt <f1> <f2> <f3> <f4>`: Apply a bandpass filter in frequency domain to the data before deconvolution. The list or tuple defines the four corner frequencies (f1, f2, f3, f4) of a cosine taper which is one between f2 and f3 and tapers to zero for f1 < f < f2 and f3 < f < f4.
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>

## 3 - Python functions
### 3.1 - Main
```python
main(args)
```
This function validates input options and decides what to do next. It also performs the parsing of `stationxml` file to get and **Inventory** object containing *network*,*station* and *channel response* informations. If the *archive mode* is activated, this function uses the `explore()` function to gather the same informations that would be parsed from a `stationxml` file. The `args` arguments contains the informations passed as options in the command line.
<br>
<br>

### 3.2 - Info
```python
help()
```
This function prints an help message with a short description of all the options available
<br>
<br>
```python
version()
```
This function prints a message with the script's version
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>


### 3.3 - Files management
```python
search_files_and_load_stream(rootdir, start_time, end_time)
```
This function opens files specified in `rootdir` and performs a **time slice** of the traces according to `start_time` and `end_time`.
<br>
<br>

```python
import_stream_from_sds_dir(rootdir, network, station, start_time, end_time)
```
This function searches for trace data starting from specified `rootdir` and associated, for **all the channels that contains waveforms**, with `network` and `station` specified in the **StationXML** file. The samples of the waveforms will be gathered according to `start_time` and `end_time`.
Please note that `rootdir` has to be a path that is compliant with SDS format: `rootdir/year/network_code/station_code/channels_dirs`
Otherwise the won't be able to correctly gather traces data.
<br>
<br>

```python
setup_output(inventory, output_path)
```
This function checks if the path specified in `output_path` exists and, if not, creates the directory and some subdirectories as follows:
```
output_path_directory
|-- deconvoluted_tracks
|-- ppsd
|   |-- npz
|   |-- plot
|-- psd
|-- raw_traces
|-- start_time.net_code.station_code.txt
```  
**NOTE:** the inventory argument was used in previous versions of the script but now is not used anymore, it will be remove in future versions.
<br>
<br>
<br>

### 3.4 - Metadata gathering

```python
explore(version, keys_list)
```
This function uses the keys_list and version informations to search recursively in the NRL Iris catalog to retrieve sensors and dataloggers responses that will be used to elaborate traces data.
<br>
<br>

### 3.5 - Elaboration and export

```python
apply_deconvolution(inv, traces, unit, prefilt)
```
This function uses informations supplied by the **Inventory object** stored in `inv` gathered from the `StationXML` file to elaborate loaded `traces` giving as result a stream with deconvoluted and converted to `unit` traces with that can be plotted or exported to an ASCII file. The `prefilt` arguments is a list of 4 values that are applied to the deconvolution as a prefilter as described above in the `--prefilt` description.
<br>
<br>
<br>

```python
export_to_ascii(output_path, start_time, end_time, traces, unit, ascii_semicolon)
```
This function exports deconvoluted traces to an ASCII file. The output will look like this:
```
Time [s]	LOCATION_CODE - CHANNEL_CODE [Export_unit];
0.00000		-0.00000008;
0.00500		-0.00000009;
0.01000		-0.00000012;
```
If the `ascii_semicolon` is set as True, the semicolon at the end of the lines and will be removed.
<br>
<br>
<br>
<br>
<br>
<br>
<br>

```python
calculate_and_plot_ppsd(inv, traces, start_t, ppsd_path, structural_mode)
```
This function uses informations supplied by the **Inventory object** stored in `inv` gathered from the `StationXML` file to elaborate loaded `traces` giving as result a `.npz` file containing traces' PPSDs that will also be plotted and stored in the specified `ppsd_path`. The `start_t` argument is used as common base through various function to give files name with the following structure:
- `start_time.TYPE_OF_EXPORT.extension` if the file contains informations about ALL the traces (**e.g.** raw traces' plots)
- `start_time.NET_CODE.STATION_CODE.LOCATION_CODE.CHANNEL_NAME.extension` if the file contains informations about a single trace (**e.g.** PPSD's plots).

If `structural_mode` arguments is set as True, parameters will be changed to produce a plot suitable for **structural analysis**, otherwise the plot will be performed set for **seismic analysis**.
<br>
<br>

```python
calculate_and_plot_psd(traces, output_unit, start_t, psd_path, export_ascii, seismic, conv_unit, ascii_semicolon)
```
This function takes a stream stored in `traces`, iterates over it, and, for each trace, simply puts together the full file path using `start_t` and `psd_path` which will be passed to the next function together with the unit label gathered from `output_unit`. If the `conv_unit` parameter is specified, il will be used to apply the conversion factor to the traces and apply the label instead of the one provided with `output_unit`.
If `seismic` arguments is set as True, the psd will be calculated using `fft_welch_seismic` instead of `fft_welch` to produce a plot suitable for **seismic analysis**.
If `export_ascii` arguments is set as True, the psd will be exported in an ASCII file with the same format of the `eport_ascii` function. The `ascii_semicolon` has the same function as described above.
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>

```python
psd_ascii(freq, filename, conv_unit, ascii_semicolon)
```
This function exports PSD values to an ASCII file. The output will look like this:
```
Time [s]	LOCATION_CODE - CHANNEL_CODE [Export_unit]
0.00000		-0.00000008
0.00500		-0.00000009
0.01000		-0.00000012
```
If the `conv_unit` parameter is specified, il will be used to apply the conversion factor to the traces and apply the label as unit of measure.
<br>
<br>

```python
fft_welch(trace, filename, unit, conv_unit)
```
This function takes a single `trace`, calculates and plots it's **PSD** using `unit` arguments as a label and stores the result in a *pdf* file according to `filename`, which is a full path.
<br>
<br>

```python
fft_welch_seismic(trace, filename, unit)
```
This function takes a single `trace`, calculates and plots it's **PSD** using `unit` arguments as a label and stores the result in a *pdf* file according to `filename`, which is a full path. If the `conv_unit` parameter is specified, il will be used to apply the conversion factor to the traces and apply the label instead of the one provided with `unit`.