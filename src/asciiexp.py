import os
from obspy import Stream
from obspy.core import UTCDateTime

unit_dictionary = {
    "ACC": "m/s^2",
    "VEL": "m/s",
    "DISP": "m"
}

output_dictionary = {
    "m/s**2": "ACC",
    "m/s": "VEL",
    "m": "DISP",
}

# CODICE LEGACY DA RIMUOVERE IN VERSIONI SUCCESSIVE
# def parallel_columns(traces, start_time, end_time, units_map, mode, ascii_semicolon, output_path, conversion_map=None, print_output=False):
    # completion = 0
    # perc_step = 100 / len(traces[0])

    # if print_output:
    #     print("Exporting to ASCII, this may take a while ... %d%%" % completion, end="", flush=True)

    # sampling_rate = traces[0].stats.sampling_rate
    # max_starttime = traces[0].stats.starttime
    # min_endtime = traces[0].stats.endtime

    # if start_time > max_starttime:
    #     max_starttime = start_time

    # if end_time < min_endtime:
    #     min_endtime = end_time

    # for tr in traces[1:]:
    #     if(tr.stats.sampling_rate != sampling_rate):
    #         if print_output:
    #             print("ERROR: Traces have different sampling rates, this feature is not supported, could not export to ASCII File")
    #         return False
    #     else:
    #         step = 100 / len(tr)
    #         if step < perc_step:
    #             perc_step = step

    #         if tr.stats.starttime > max_starttime:
    #             max_starttime = tr.stats.starttime
    #         if tr.stats.endtime < max_starttime:
    #             min_endtime = tr.stats.endtime

    # time_str = max_starttime.strftime("%Y%m%d%H%M%SS")

    # file_name = os.path.join(output_path, time_str + ".ASCII.txt")

    # f = open(file_name, "w+")

    # header_line1 = "\t\tASCII Export - Lunitek\n"
    # f.write(header_line1)

    # header_line2 = "Time[s]"
    # for tr in traces:

    #     if mode == "NATURAL":
    #         for ch_map in units_map:
    #             if tr.stats.location == ch_map['loc'] and tr.stats.channel == ch_map['channel']:
    #                 calc_output = output_dictionary[ch_map['unit'].lower()]
    #                 unit = ch_map['unit']
    #                 break
    #     else:
    #         calc_output = mode
    #         unit = unit_dictionary[calc_output]

    #     if conversion_map:
    #         for conv in conversion_map:
    #             if calc_output == conv["source"]:
    #                 unit = conv["label"]

    #     header_line2 += "\t"
    #     header_line2 += tr.stats.network
    #     header_line2 += " - "
    #     header_line2 += tr.stats.station
    #     header_line2 += " - "
    #     if tr.stats.location:
    #         header_line2 += tr.stats.location
    #         header_line2 += " - "
    #     header_line2 += tr.stats.channel
    #     header_line2 += "["
    #     header_line2 += unit
    #     header_line2 += "]"
    # if not ascii_semicolon:
    #     header_line2 += ";"
    # header_line2 += "\n"

    # f.write(header_line2)

    # current_line_time = max_starttime

    # i = 0
    # while current_line_time < min_endtime:
    #     time = float(current_line_time - max_starttime)
    #     line = '%.5f' % time
    #     for tr in traces:
    #         line += "\t"
    #         try:
    #             line += '%.8f' % tr.data[i]
    #         except IndexError:
    #             line += '%.8f' % 0.00000000
    #     if not ascii_semicolon:
    #         line += ";"
    #     f.write(line + "\n")

    #     current_line_time += (1 / sampling_rate)
    #     i += 1
    #     completion += perc_step
    #     if print_output:
    #         print("\rExporting to ASCII, this may take a while ... %d%%" % completion, end="", flush=True)
    # if print_output:
    #     print("\rExporting to ASCII, this may take a while ... Done")
    # return True


def time_export(
    traces: Stream,
    start_time: UTCDateTime,
    end_time: UTCDateTime,
    units_map: list,
    mode: str,
    output_path: str,
    print_header: bool,
    separator: str,
    endline: str,
    numformat: dict,
    print_xaxis: bool,
    conversion_map: list = None,
    pretrigger: int = 0,
    print_output: bool = False
) -> bool:

    completion = 0
    perc_step = 100 / len(traces[0])

    if print_output:
        print("Exporting to ASCII, this may take a while ... %d%%" % completion, end="", flush=True)

    sampling_rate = traces[0].stats.sampling_rate
    max_starttime = traces[0].stats.starttime
    min_endtime = traces[0].stats.endtime

    if start_time > max_starttime:
        max_starttime = start_time

    if end_time < min_endtime:
        min_endtime = end_time

    for tr in traces[1:]:
        if(tr.stats.sampling_rate != sampling_rate):
            if print_output:
                print("ERROR: Traces have different sampling rates, this feature is not supported, could not export to ASCII File")
            return False
        else:
            step = 100 / len(tr)
            if step < perc_step:
                perc_step = step

            if tr.stats.starttime > max_starttime:
                max_starttime = tr.stats.starttime
            if tr.stats.endtime < max_starttime:
                min_endtime = tr.stats.endtime

    time_str = max_starttime.strftime("%Y%m%d%H%M%SS")

    file_name = os.path.join(output_path, time_str + ".ASCII.txt")

    f = open(file_name, "w+")

    header_line1 = "\t\tASCII Export - Lunitek\n"

    if print_header:
        f.write(header_line1)

    header_line2 = ""
    if print_xaxis:
        header_line2 += "Time[s]"
        # header_line2 += "\t"
        header_line2 += separator

    for k, tr in enumerate(traces):
        if mode == "NATURAL":
            for ch_map in units_map:
                if tr.stats.location == ch_map['loc'] and tr.stats.channel == ch_map['channel']:
                    calc_output = output_dictionary[ch_map['unit'].lower()]
                    unit = ch_map['unit']
                    break
        else:
            calc_output = mode
            unit = unit_dictionary[calc_output]

        if conversion_map:
            for conv in conversion_map:
                if calc_output == conv["source"]:
                    unit = conv["label"]

        header_line2 += tr.stats.network
        header_line2 += " - "
        header_line2 += tr.stats.station
        header_line2 += " - "
        if tr.stats.location:
            header_line2 += tr.stats.location
            header_line2 += " - "
        header_line2 += tr.stats.channel
        header_line2 += "["
        header_line2 += unit
        header_line2 += "]"
        if k != len(traces) - 1:
            header_line2 += separator

    header_line2 += endline
    # header_line2 += ";"
    # header_line2 += "\n"

    if print_header:
        f.write(header_line2)

    current_line_time = max_starttime - pretrigger

    i = 0
    while current_line_time < (min_endtime - pretrigger):
        time = float(current_line_time - max_starttime)
        line = ""
        if print_xaxis:
            line += numformat["x"].format(time)
            line += separator
        for j, tr in enumerate(traces):
            try:
                # line += '%.8f' % tr.data[i]
                line += numformat["y"].format(tr.data[i])
            except IndexError:
                # line += '%.8f' % 0.00000000
                line += numformat["y"].format(0.00000000)
            if j != len(traces) - 1:
                line += separator

        line += endline
        f.write(line)

        current_line_time += (1 / sampling_rate)
        i += 1
        completion += perc_step
        if print_output:
            print("\rExporting to ASCII, this may take a while ... %d%%" % completion, end="", flush=True)
    if print_output:
        print("\rExporting to ASCII, this may take a while ... Done")
    return True


def psd_old(psds, start_t, ascii_semicolon, path, print_output=False):

    completion = 0
    perc_step = 100 / len(psds)

    if print_output:
        print("Exporting PSD to ASCII, this may take a while ... %d%%" % completion, end="", flush=True)

    # filename = "" # os.path.join(output_path, time_str + ".ASCII.txt")

    for psd in psds:
        psd_file_name = start_t.strftime("%Y%m%d%H%M%S") + '.'
        psd_file_name += psd["trace"].stats.network + '.' + psd["trace"].stats.station + '.' + psd["trace"].stats.location + "." + psd["trace"].stats.channel + ".psd"
        filename = path + psd_file_name + ".txt"

        f = open(filename, "w+")

        header_line1 = "\t\tPSD ASCII Export - Lunitek\n"
        f.write(header_line1)

        header_line2 = "Frequencies[Hz]\tPSD [(" + psd["unit"] + ")^2/Hz]"

        if not ascii_semicolon:
            header_line2 += ";"
        header_line2 += "\n"

        f.write(header_line2)

        for i in range(len(psd["f"])):
            freq_str = str('%.5f' % psd["f"][i])
            psd_str = str(psd["Pxx_den"][i])
            line = '%s' % freq_str
            # if len(freq_str) < 6:
            #     line += "\t"
            line += "\t"
            line += '%s' % psd_str
            if not ascii_semicolon:
                line += ";"
            f.write(line + "\n")

        completion += perc_step

        if print_output:
            print("\rExporting PSD to ASCII, this may take a while ... %d%%" % completion, end="", flush=True)

    if print_output:
        print("\rExporting PSD to ASCII, this may take a while ... Done")


def psd(
    psds: list,
    start_t: UTCDateTime,
    print_header: bool,
    separator: str,
    endline: str,
    numformat: dict,
    print_xaxis: bool,
    output_path: str,
    print_output: bool = False
) -> None:

    completion = 0
    perc_step = 100 / len(psds)

    if print_output:
        print("Exporting PSD to ASCII, this may take a while ... %d%%" % completion, end="", flush=True)

    for psd in psds:
        psd_file_name = start_t.strftime("%Y%m%d%H%M%S") + '.'
        psd_file_name += psd["trace"].stats.network + '.' + psd["trace"].stats.station + '.' + psd["trace"].stats.location + "." + psd["trace"].stats.channel + ".psd"
        filename = output_path + psd_file_name + ".txt"

        f = open(filename, "w+")

        if print_header:
            header_line1 = "\t\tPSD ASCII Export - Lunitek\n"
            f.write(header_line1)

            header_line2 = "Frequencies[Hz]\tPSD [(" + psd["unit"] + ")^2/Hz]"
            header_line2 += endline

            f.write(header_line2)

        for i in range(len(psd["f"])):
            freq_str = numformat["x"].format(psd["f"][i])
            psd_str = numformat["y"].format(psd["Pxx_den"][i])
            line = ""
            if print_xaxis:
                line += freq_str
                line += separator
            # if len(freq_str) < 6:
            #     line += "\t"
            line += psd_str
            line += endline
            f.write(line)

        completion += perc_step

        if print_output:
            print("\rExporting PSD to ASCII, this may take a while ... %d%%" % completion, end="", flush=True)

    if print_output:
        print("\rExporting PSD to ASCII, this may take a while ... Done")
