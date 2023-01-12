import os
import obspy
import argparse
from obspy.core import UTCDateTime
from obspy.core import Stream

import nrlutils
import files
import plots
import calc
import asciiexp
import preprocessing

VERSION = "1.0"

DEFAULT_EXPORT_MODE = "NATURAL"
DEFAULT_XML_INPUT = ""
DEFAULT_EVENT_CODE = "00000000"

DEFAULT_START_TIME = "2000-01-01T00:00:00"
DEFAULT_END_TIME = "2100-01-01T00:00:00"
DEFAULT_PREFILT = [0.1, 0.2, 40, 45]
DEFAULT_PRETRIGGER = 0

DEFAULT_IRIS_PATH = "IRIS"

output_dictionary = {
    "m/s**2": "ACC",
    "m/s": "VEL",
    "m": "DISP",
    "PA": "INFRA"
}

unit_dictionary = {
    "ACC": "m/s^2",
    "VEL": "m/s",
    "DISP": "m",
    "INFRA": "PA",
    "NATURAL": ""
}



default_export_settings = {
    "ascii": [
        {
            "type": "time",
            "x_axis": True,
            "format": {
                "x": "{:.5f}",
                "y": "{:.5f}"
            },
            "header": True,
            "separator": "\t",
            "endline": ";\n"
        },
        {
            "type": "psd",
            "x_axis": True,
            "format": {
                "x": "{:.5f}",
                "y": "{:.8e}"
            },
            "header": True,
            "separator": "\t",
            "endline": ";\n\r"
        }
    ]
}


def version():
    print("LUNITEK - MSEED PLOT TOOL V" + VERSION)


def help():
    print("")
    print("\t\tLUNITEK - MSEED PLOT TOOL V%s\n" % VERSION)
    print("\t\t\tHELP GUIDE\n")
    print("--input=<path>\t\t\t\tSpecifies path where miniSEED are placed")
    print("--output=<path>\t\t\t\tSpecifies path where output files will be saved")
    print("--export=<ACC|VEL|DISP>\t\t\tSpecifies output unit for deconvolution operation")
    print("--stationxml=<file.xml>\t\t\tSpecifies StationXML file for nominal response calculation")
    print("--archive=<event_CODE.tar.gz>\t\t(Optional) Specifies tar.gz archive to work in archive mode, input path will be ignored")
    print("--start_time=<2000-01-01T00:00:00>\t(Optional) Specifies start time for time window")
    print("--end_time=<2100-01-01T00:00:00>\t(Optional) Specifies end time for time window")
    print("--plot\t\t\t\t\tDeconvoluted traces will be plotted to a PDF file in output path")
    print("--plot_ppsd\t\t\t\tPPSD of traces will be calculated, plotted to PDF files and saved in binary files in output path")
    print("--plot_psd\t\t\t\tPSD of traces will be calculated and plotted to a PDF file in output path")
    print("--plot_raw\t\t\t\tRaw traces will be plotted to a PDF file in output path")
    print("--ascii\t\t\t\t\tDeconvoluted traces will be exported to an ASCII file in output path")
    print("--psd_ascii\t\t\t\tPSD data will be exported to an ASCII file in output path")
    print("--sds\t\t\t\t\tInput path will be treated as an SDS directory")
    print("--structural\t\t\t\tPPSD Plot will be enabled and performed in Structural mode")
    print("--seismic\t\t\t\tPSD Plot will be enabled and performed in Seismic mode")
    print("--infrasound\t\t\t\tPSD Plot will be enabled and performed in Infrasound mode")
    print("--spectrogram\t\t\t\tThe spectrogram of the input stream will be calculated and plotted")
    print("--dayplot\t\t\t\tA day plot of the traces will be performed")
    print("--print_output\t\t\t\tThe script will print on stdout current progress")
    print("--iris=<path>\t\t\t\tSets IRIS directory path")
    print("--unit <source> <label> <factor> <source> <label> <factor> <source> <label> <factor>\tTraces will be converted by multiplying values for <factor> value and some plots/exports will use converted traces")
    print("--ascii_semicolon\t\t\tThe semicolon won't be used as end-of-line character in ASCII eport files")
    print("--prefilt <f1> <f2> <f3> <f4>\t\tApply a bandpass filter in frequency domain to the data before deconvolution.\n"
    "\t\t\t\t\tThe list or tuple defines the four corner frequencies (f1, f2, f3, f4) of a cosine taper\n"
    "\t\t\t\t\twhich is one between f2 and f3 and tapers to zero for f1 < f < f2 and f3 < f < f4.")
    print("--channels_list\t\t\t\tSpecify explicit channels list to be converted")
    print("--settings\t\t\t\tFile with conversion and output settings")
    print("--fill_missing\t\t\t\tFills ASCII file with zeros if some specified channels are missing in input data")
    print("--reorder\t\t\t\tReorders data stream according to specified channels list")
    print("--flat_convert\t\t\t\tConverts data using channels sensitivity instead of applying deconvolution")
    print("--rms\t\t\t\t\tCalculates RMS on detrended channels data")
    print("-h,--help\t\t\t\tPrints this help")
    print("-v,--version\t\t\t\tPrints script's version")
    print("")


def main(args):
    input_path = ""
    stationxml = ""
    nrl_catalog = nrlutils.get_nrl_catalog(args.iris)
    channels_map = []

    conversion_map = []

    if not args.archive and not args.input_path:
        print("ERROR: Input path not specified")
        exit(-1)
    elif args.input_path:
        input_path = args.input_path

    # VERIFY TIME WINDOW
    try:
        start_t = UTCDateTime(args.start_t_str)
        end_t = UTCDateTime(args.end_t_str)
    except TypeError:
        print("ERROR: Cannot create time window")
        exit(-1)

    if start_t > end_t:
        if args.print_output:
            print("WARNING: Start time can't be greater than End time, I'm swapping them")
        tmp_t = start_t
        start_t = end_t
        end_t = tmp_t

    # LOAD CONVERSION MAP
    if args.conversion:
        if len(args.conversion) % 3 != 0 or len(args.conversion) > 9:
            print("ERROR: Not a valid conversion list")
            exit(-1)
        try:
            factor = float(args.conversion[2])
        except Exception as e:
            print(e)
            print("ERROR: Not a valid conversion factor")
            exit(-1)

        first_map = {
            "source": args.conversion[0],
            "label": args.conversion[1],
            "factor": factor
        }
        conversion_map.append(first_map)
        if len(args.conversion) > 3:
            try:
                factor = float(args.conversion[5])
            except Exception as e:
                print(e)
                print("ERROR: Not a valid conversion factor")
                exit(-1)

            second_map = {
                "source": args.conversion[3],
                "label": args.conversion[4],
                "factor": factor
            }
            conversion_map.append(second_map)
        if len(args.conversion) > 6:
            try:
                factor = float(args.conversion[8])
            except Exception as e:
                print(e)
                print("ERROR: Not a valid conversion factor")
                exit(-1)

            third_map = {
                "source": args.conversion[6],
                "label": args.conversion[7],
                "factor": factor
            }
            conversion_map.append(third_map)

    # READ PRETRIGGER
    try:
        pretrigger = int(args.pretrigger)
    except Exception:
        print("ERROR: Could not parse pretrigger")
        exit(-1)

    # LOAD STATIONXML INFORMATIONS
    if args.archive:
        if args.print_output:
            print("INFO: Working in archive mode, extracting data and fetching StationXML, please wait just a moment..")
            print("")
        input_path, poseidon_config_file, factory_data_file = files.init_archive_mode(args.archive)
        stationxml, channels_map, pretrigger = nrlutils.stationxml_from_metadata(poseidon_config_file, factory_data_file, nrl_catalog)
    else:
        if not args.stationxml:
            print("ERROR: StationXML file not specified")
            exit(-1)
        if not os.path.isfile(args.stationxml):
            print("ERROR: Cannot find StationXML file")
            exit(-1)
        else:
            stationxml = args.stationxml

    # DELETE PRETRIGGER IF NOT WANTED
    if not args.print_pretrigger:
        pretrigger = 0

    inv = obspy.read_inventory(stationxml)

    # UPDATE CHANNELS MAP WITH UNITS
    for inv_ch in inv[0][0]:
        unit_map = {}
        unit_map['loc'] = inv_ch.location_code
        unit_map['ch'] = inv_ch.code

        for i, ch in enumerate(channels_map):
            if ch["loc"] == unit_map["loc"] and ch["channel"] == unit_map["ch"]:
                channels_map[i]["unit"] = inv_ch.response.response_stages[0].input_units

    if args.print_output:
        print("")
        print("\t\tRECAP")
        print("START TIME:\t\t\t%s" % start_t.strftime("%Y-%m-%dT%H-%M-%S"))
        print("END TIME:\t\t\t%s" % end_t.strftime("%Y-%m-%dT%H-%M-%S"))
        print("STATION XML FILE:\t\t%s" % stationxml)
        print("IRIS PATH:\t\t\t%s" % args.iris)
        if args.conversion:
            print("UNIT CONVERSION:\t\t%s\t%s" % (args.conversion[0], args.conversion[1]))

    # SETUP OUTPUT DIRECTORY
    if not args.output_path:
        output_base_path = "Exports"
    else:
        output_base_path = args.output_path
    if args.print_output:
        print("INPUT PATH:\t\t\t%s" % input_path)
        print("OUTPUT PATH:\t\t\t%s" % output_base_path)

    if not output_base_path.endswith("/"):
        output_base_path += "/"

    err = files.setup_output(output_base_path)
    if type(err) is tuple:
        print(err[1])
        exit(err[0])

    raw_traces_path = output_base_path + "raw_traces/"
    elab_traces_path = output_base_path + "deconvoluted_traces/"
    ppsd_path = output_base_path + "ppsd/"
    psd_path = output_base_path + "psd/"
    ascii_path = output_base_path + "ASCII/"
    spectrogram_path = output_base_path + "spectrogram/"
    dayplot_path = output_base_path + "dayplot/"

    # READ OUTPUT UNIT
    export_mode = args.export_mode

    if export_mode not in unit_dictionary:
        print("ERROR: Invalid export type %s" % export_mode)
        exit(-1)
    if args.print_output:
        print("EXPORT TYPE:\t\t\t%s" % export_mode)

        # PLOT MESSAGES
        if args.plot_traces:
            print("PLOT DECONVOLUTED TRACES:\tEnabled")
        else:
            print("PLOT DECONVOLUTED TRACES:\tDisabled")

        if args.plot_raw_traces:
            print("PLOT RAW TRACES:\t\tEnabled")
        else:
            print("PLOT RAW TRACES:\t\tDisabled")
        if args.ascii_exp:
            print("ASCII EXPORT:\t\t\tEnabled")
        else:
            print("ASCII EXPORT:\t\t\tDisabled")

        if args.plot_psd:
            print("PLOT PSD:\t\t\tEnabled")
        else:
            print("PLOT PSD:\t\t\tDisabled")

        if args.plot_ppsd or args.structural_mode:
            print("PLOT PPSD:\t\t\tEnabled")
        else:
            print("PLOT PPSD:\t\t\tDisabled")

        if args.plot_ppsd and not args.structural_mode:
            print("PPSD MODE:\t\t\tSeismic")
        elif args.structural_mode:
            print("PPSD MODE:\t\t\tStructrural")

        if args.plot_spectrogram:
            print("PLOT SPECTROGRAM:\t\tEnabled")
        else:
            print("PLOT SPECTROGRAM:\t\tDisabled")

        # LOAD DATA STREAMS
        print("\n\tTRACES INFO\n")

    if args.sds_search:
        loaded_stream = Stream()
        for net in inv.networks:
            for sta in net.stations:
                loaded_stream += files.import_stream_from_sds_dir(input_path, net.code, sta.code, start_t, end_t, args.print_output)
    else:
        if os.path.isdir(input_path):
            if not input_path.endswith("/"):
                input_path += "/"
            input_path += "*"
        loaded_stream = files.search_files_and_load_stream(input_path, start_t, end_t, args.print_output)

    if type(loaded_stream) is tuple:
        print(loaded_stream[1])
        exit(loaded_stream[0])

    if len(loaded_stream) < 1:
        print("ERROR: No traces found")
        return -1

    start_t_tmp = loaded_stream[0].stats.starttime
    end_t_tmp = loaded_stream[0].stats.endtime

    for tr in loaded_stream[1:]:
        if tr.stats.starttime > start_t_tmp:
            start_t_tmp = tr.stats.starttime
        if tr.stats.endtime < end_t:
            end_t_tmp = tr.stats.endtime

    if start_t < start_t_tmp:
        start_t = start_t_tmp

    if end_t > end_t_tmp:
        end_t = end_t_tmp

    # START ELABORATION
    if args.print_output:
        print("\tELABORATION STARTED\n")

    # RAW TRACES
    traces = loaded_stream.copy()

    # MIN - MAX
    # traces.merge()
    # import numpy as np
    # for trace in traces:
    #     name = trace.stats.network + "." + trace.stats.station + "." + trace.stats.location + "." + trace.stats.channel
    #     maxval = trace.data.max()
    #     minval = trace.data.min()
    #     print("{}\tMAX: {}\tMIN: {}".format(name, maxval, minval))

    # LOAD EXPORT SETTINGS

    if args.export_settings:
        export_settings = files.load_settings_file(args.export_settings)
    else:
        export_settings = default_export_settings

    ERROR = False

    # PLOT RAW TRACES
    if args.plot_raw_traces:
        plots.raw(traces, start_t, raw_traces_path, args.print_output)

    # PLOT DAYPLOT
    if args.dayplot:
        plots.dayplot(traces, start_t, dayplot_path, args.print_output)

    for i, elem in enumerate(args.prefilt):
        args.prefilt[i] = float(elem)

    if args.calc_rms:
        calc.detrended_rms(traces, inv)


    if args.flat_convert:
        elab_traces = Stream()
        for trace in traces:
            for i, network in enumerate(inv):
                if network.code == trace.stats.network:
                    net = i
                    for j, station in enumerate(inv[net]):
                        if station.code == trace.stats.station:
                            sta = j
                            for k, channel in enumerate(inv[net][sta]):
                                if channel.code == trace.stats.channel and channel.location_code == trace.stats.location:
                                    sensitivity = channel.response.instrument_sensitivity.value
                                    print(sensitivity)
                                    elab_trace = trace.copy()
                                    elab_trace.data = (elab_trace.data / sensitivity)
                                    elab_traces.append(elab_trace)
    # REMOVE RESPONSE
    else:
        elab_traces = calc.apply_deconvolution(inv, traces, channels_map, export_mode, args.prefilt, args.print_output)  # DA CAMBIARE

    # CONVERSION
    if args.conversion:
        conv_traces = []
        calc_output = args.export_mode
        for trace in elab_traces:
            for ch_map in channels_map:
                if trace.stats.location == ch_map['loc'] and trace.stats.channel == ch_map['channel']:
                    calc_output = output_dictionary[ch_map['unit'].lower()]
                    break
            for conv in conversion_map:
                if calc_output == conv["source"]:
                    # trace.data *= conv["factor"]
                    trace.data = (trace.data * conv["factor"])

    if args.plot_traces:
        plots.deconvoluted(elab_traces, channels_map, export_mode, start_t, elab_traces_path, args.print_output)

    # CALCULATE AND PLOT PPSD
    if args.plot_ppsd or args.structural_mode:
        ppsds = calc.ppsd(inv, loaded_stream, start_t, args.structural_mode, args.print_output)
        ERROR = plots.ppsd(ppsds, args.structural_mode, ppsd_path, args.print_output)

    if args.plot_psd or args.psd_ascii or args.seismic_mode:
        if args.conversion:
            psds = calc.psd(elab_traces, channels_map, export_mode, conversion_map)
        else:
            psds = calc.psd(elab_traces, channels_map, export_mode)

        if args.seismic_mode:
            ERROR = plots.psd_seismic(psds, channels_map, export_mode, start_t, psd_path)
        elif args.infrasound_mode:
            ERROR = plots.psd_infrasound(psds, channels_map, export_mode, start_t, psd_path)
        else:
            ERROR = plots.psd(psds, channels_map, export_mode, start_t, psd_path)

        if args.psd_ascii:
            psd_idx = files.get_index_with_key_value(export_settings['ascii'], "type", "psd")
            psd_print_header = export_settings['ascii'][psd_idx]['header']
            psd_separator = export_settings['ascii'][psd_idx]['separator']
            psd_endline = export_settings['ascii'][psd_idx]['endline']
            psd_numformat = export_settings['ascii'][psd_idx]['format']
            psd_print_xaxis = export_settings['ascii'][psd_idx]['x_axis']

            # asciiexp.psd(psds, start_t, args.ascii_semicolon, psd_path + "ASCII/", args.print_output)
            asciiexp.psd(psds, start_t, psd_print_header, psd_separator, psd_endline, psd_numformat, psd_print_xaxis, psd_path + "ASCII/", args.print_output)

    # PLOT SPECTROGRAM
    if args.plot_spectrogram:
        sp_traces = elab_traces
        if args.unit:
            sp_traces = conv_traces
        plots.spectrogram(sp_traces, spectrogram_path, args.print_output)

    # ASCII EXPORT
    if args.ascii_exp:
        # ERROR = not asciiexp.parallel_columns(elab_traces, start_t, end_t, channels_map, export_mode, args.ascii_semicolon, ascii_path, conversion_map, args.print_output) or ERROR
        time_idx = files.get_index_with_key_value(export_settings['ascii'], "type", "time")
        time_print_header = export_settings['ascii'][time_idx]['header']
        time_separator = export_settings['ascii'][time_idx]['separator']
        time_endline = export_settings['ascii'][time_idx]['endline']
        time_numformat = export_settings['ascii'][time_idx]['format']
        time_print_xaxis = export_settings['ascii'][time_idx]['x_axis']

        channels_list = []

        if "traces" in export_settings['ascii'][time_idx]:
            channels_list = export_settings['ascii'][time_idx]['traces']

        if args.channels_list:
            channels_list = args.channels_list

        if args.fill_missing:
            if not channels_list:
                print("ERROR: Missing channels' list, please specify channels")
                exit(-1)
            preprocessing.fill_stream(elab_traces, channels_list)
        if args.reorder:
            if not channels_list:
                print("ERROR: Missing channels' list, please specify channels")
                exit(-1)
            elab_traces = preprocessing.reorder_stream(elab_traces, channels_list)

        ERROR = not asciiexp.time_export(elab_traces, start_t, end_t, channels_map, export_mode, ascii_path, time_print_header, time_separator, time_endline, time_numformat, time_print_xaxis, conversion_map, pretrigger, args.print_output) or ERROR

    if args.print_output:
        print("")
        if not ERROR:
            print("All tasks completed successfully, exiting now, bye bye!")
        else:
            print("Not all tasks could be completed, please check error messages")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="mSEED data processing", add_help=False)

    parser.add_argument("--input", action="store", dest="input_path")
    parser.add_argument("--output", action="store", dest="output_path")
    parser.add_argument("--start_time", action="store", dest="start_t_str", default=DEFAULT_START_TIME)
    parser.add_argument("--end_time", action="store", dest="end_t_str", default=DEFAULT_END_TIME)
    parser.add_argument("--export", action="store", dest="export_mode", default=DEFAULT_EXPORT_MODE)
    parser.add_argument("--stationxml", action="store", dest="stationxml", default=DEFAULT_XML_INPUT)
    parser.add_argument("--archive", action="store", dest="archive")
    parser.add_argument("--plot", action="store_true", dest="plot_traces", default=False)
    parser.add_argument("--plot_raw", action="store_true", dest="plot_raw_traces", default=False)
    parser.add_argument("--plot_ppsd", action="store_true", dest="plot_ppsd", default=False)
    parser.add_argument("--plot_psd", action="store_true", dest="plot_psd", default=False)
    parser.add_argument("--ascii", action="store_true", dest="ascii_exp", default=False)
    parser.add_argument("--psd_ascii", action="store_true", dest="psd_ascii", default=False)
    parser.add_argument("--sds", action="store_true", dest="sds_search", default=False)
    parser.add_argument("--structural", action="store_true", dest="structural_mode", default=False)
    parser.add_argument("--seismic", action="store_true", dest="seismic_mode", default=False)
    parser.add_argument("--infrasound", action="store_true", dest="infrasound_mode", default=False)
    parser.add_argument("--spectrogram", action="store_true", dest="plot_spectrogram", default=False)
    parser.add_argument("--dayplot", action="store_true", dest="dayplot", default=False)
    parser.add_argument("--print_output", action="store_true", dest="print_output", default=False)
    parser.add_argument("--iris", action="store", dest="iris", default=DEFAULT_IRIS_PATH)
    parser.add_argument("--unit", action="store", dest="conversion", nargs='+')
    parser.add_argument("--ascii_semicolon", action="store_true", dest="ascii_semicolon", default=False)
    parser.add_argument("--prefilt", action="store", dest="prefilt", nargs=4, default=DEFAULT_PREFILT)
    parser.add_argument("--fill_missing", action="store_true", dest="fill_missing", default=False)
    parser.add_argument("--reorder", action="store_true", dest="reorder", default=False)
    parser.add_argument("--channels_list", action="store", dest="channels_list", nargs='+', default=False)
    parser.add_argument("--settings", action="store", dest="export_settings")
    parser.add_argument("--flat_convert", action="store_true", dest="flat_convert", default=False)
    parser.add_argument("--rms", action="store_true", dest="calc_rms", default=False)
    parser.add_argument("--print_pretrigger", action="store_true", dest="print_pretrigger", default=False)
    parser.add_argument("--pretrigger", action="store", dest="pretrigger", default=DEFAULT_PRETRIGGER)
    parser.add_argument("-h", "--help", action="store_true", dest="help", default=False)
    parser.add_argument("-v", "--version", action="store_true", dest="version", default=False)
    args = parser.parse_args()

    args.print_output = args.print_output

    if args.help:
        help()
        exit(0)
    if args.version:
        version()
        exit(0)
    else:
        main(args)
