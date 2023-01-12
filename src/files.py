import shutil
import os
import obspy
from obspy.clients.filesystem.sds import Client as SDSClient
import json


def init_archive_mode(archive_file):
    dirname_parts = archive_file.split("/")
    dirname = dirname_parts[-1]  # .split(".")[0]
    dirname = dirname.replace(".tar.gz", "")
    shutil.unpack_archive(archive_file, dirname)
    os.chdir(dirname)

    mseed_folder = "mseed"

    try:
        os.stat(mseed_folder)
        shutil.rmtree(mseed_folder)
        os.mkdir(mseed_folder)
    except Exception:
        os.mkdir(mseed_folder)

    mseed_files = os.listdir(".")

    for f in mseed_files:
        if f.endswith(mseed_folder):
            shutil.move(f, mseed_folder)

    poseidon_config = ""
    factory_data = ""

    if os.path.isfile("poseidon_config.xml"):
        poseidon_config = os.path.abspath("poseidon_config.xml")
    else:
        print("ERROR: poseidon_config.xml not present in specified archive")
        exit(-1)

    if os.path.isfile("factory_data.xml"):
        factory_data = os.path.abspath("factory_data.xml")
    else:
        print("ERROR: factory_data.xml not present in specified archive")
        exit(-1)

    return mseed_folder, poseidon_config, factory_data


def setup_output(output_path):
    try:
        if not os.path.exists(output_path):
            os.mkdir(output_path)
            os.mkdir(output_path + "deconvoluted_traces")
            os.mkdir(output_path + "ASCII")
            os.mkdir(output_path + "raw_traces")
            os.mkdir(output_path + "psd")
            os.mkdir(output_path + "psd/ASCII")
            os.mkdir(output_path + "ppsd")
            os.mkdir(output_path + "ppsd/npz")
            os.mkdir(output_path + "ppsd/plot")
            os.mkdir(output_path + "spectrogram")
            os.mkdir(output_path + "dayplot")
        else:
            if not os.path.exists(output_path + "deconvoluted_traces"):
                os.mkdir(output_path + "deconvoluted_traces")
            if not os.path.exists(output_path + "raw_traces"):
                os.mkdir(output_path + "raw_traces")
            if not os.path.exists(output_path + "ASCII"):
                os.mkdir(output_path + "ASCII")
            if not os.path.exists(output_path + "psd"):
                os.mkdir(output_path + "psd")
            if not os.path.exists(output_path + "psd/ASCII"):
                os.mkdir(output_path + "psd/ASCII")
            if not os.path.exists(output_path + "ppsd"):
                os.mkdir(output_path + "ppsd")
            if not os.path.exists(output_path + "ppsd/npz"):
                os.mkdir(output_path + "ppsd/npz")
            if not os.path.exists(output_path + "ppsd/plot"):
                os.mkdir(output_path + "ppsd/plot")
            if not os.path.exists(output_path + "spectrogram"):
                os.mkdir(output_path + "spectrogram")
            if not os.path.exists(output_path + "dayplot"):
                os.mkdir(output_path + "dayplot")

        return 0
    except OSError:
            errmsg = "ERROR: Cannot create output path"
            return (-1, errmsg)


def import_stream_from_sds_dir(rootdir, network, station, start_time, end_time, print_output=False):
    client = SDSClient(rootdir)
    stream = client.get_waveforms(network, station, "*", "*", start_time, end_time)
    if print_output:
        print("Traces found:")
        for tr in stream:
            print(tr.stats)
            print("")
    return stream


def search_files_and_load_stream(rootdir, start_time, end_time, print_output=False):
    try:
        tmp_stream = obspy.read(rootdir)
    except Exception:
        errmsg = "ERROR: Cannot read traces from specified input path. Please check file presence and make sure that it is a leaf path (no directories in it)"
        return (-1, errmsg)

    stream = obspy.Stream()

    deny_channel_list = ['SOH', 'LOG', 'RNG']

    for tr in tmp_stream:
        if tr.stats.channel not in deny_channel_list:
            stream += tr

    # print("LUNGHEZZA STREAM %d" % len(stream))

    stream.merge(method=1)

    max_starttime = stream[0].stats.starttime
    min_endtime = stream[0].stats.endtime

    for tr in stream[1:]:
        if tr.stats.starttime > max_starttime:
            max_starttime = tr.stats.starttime
        if tr.stats.endtime < max_starttime:
            min_endtime = tr.stats.endtime

    if start_time > max_starttime and start_time < min_endtime:
        max_starttime = start_time

    if end_time < min_endtime and end_time > max_starttime:
        min_endtime = end_time

    # print(max_starttime)
    # print(min_endtime)

    for i in range(len(stream)):
        stream[i] = stream[i].slice(max_starttime, min_endtime)

    if print_output:
        print("Traces found:")
        for tr in stream:
            print(tr.stats)
            print("")

    return stream


def load_settings_file(path):
    abspath = os.path.abspath(path)
    with open(abspath) as settings_file:
        export_settings = json.load(settings_file)
        return export_settings


def get_index_with_key_value(dict_list, key, value):
    for i, elem in enumerate(dict_list):
        if key in elem and elem[key] == value:
            return i
    return None
