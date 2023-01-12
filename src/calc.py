from obspy.core import Stream
from obspy.signal import PPSD

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

# [0.1, 0.2, 40, 45]


def apply_deconvolution(inv, traces, units_map, mode, prefilt=[0, 0, 0, 0], print_output=False):

    completion = 0
    perc_step = 100 / len(traces)
    if print_output:
        print("Applying deconvolution, this may take a while ... %d%%" % completion, end="", flush=True)
    for trace in traces:
        if mode == "NATURAL":
            for ch_map in units_map:
                if trace.stats.location == ch_map['loc'] and trace.stats.channel == ch_map['channel']:
                    calc_output = output_dictionary[ch_map['unit'].lower()]
                    break
        else:
            calc_output = mode
        trace.remove_response(inventory=inv, water_level=1000, pre_filt=None, output=calc_output)
        completion += perc_step
        if print_output:
            print("\rApplying deconvolution, this may take a while ... %d%%" % completion, end="", flush=True)
    if print_output:
        print("\rApplying deconvolution, this may take a while ... Done")
    return traces


def ppsd(inv, traces, start_t, structural_mode, print_output=False):
    completion = 0
    perc_step = 100 / len(traces)

    ppsds = []
    if print_output:
        print("Calculating and plotting PPSDs, this may take a while ... %d%%" % completion, end="", flush=True)

    for trace in traces:
        st = Stream(trace)
        ppsd_file_name = start_t.strftime("%Y%m%d%H%M%S") + '.' + trace.stats.network + '.' + trace.stats.station + '.' + trace.stats.location + "." + trace.stats.channel + ".ppsd"
        if structural_mode:
            ppsd = PPSD(st[0].stats, metadata=inv, ppsd_lenght=10.0)
        else:
            ppsd = PPSD(st[0].stats, metadata=inv)
        ppsd.add(st)
        new_ppsd = {
            "name": ppsd_file_name,
            "data": ppsd
        }
        ppsds.append(new_ppsd)

        completion += perc_step
        if print_output:
            print("\rCalculating and plotting PPSDs, this may take a while ... %d%%" % completion, end="", flush=True)
    if print_output:
        print("\rCalculating and plotting PPSDs, this may take a while ... Done")

    return ppsds


def fft_welch(trace, conversion_factor, perseg=20000, overlap=10000):
    from scipy import signal

    x = trace.data
    fs = trace.stats.sampling_rate

    nperseg = perseg
    noverlap = overlap
    if nperseg > len(x):
        nperseg = len(x)
        noverlap = nperseg / 2

    f, Pxx_den = signal.welch(x, fs, nperseg=nperseg, noverlap=noverlap)

    Pxx_den = Pxx_den * conversion_factor

    return f, Pxx_den


def psd(traces, units_map, export_mode, conversion_map=None):
    psds = []

    for trace in traces:
        if export_mode == "NATURAL":
            for ch_map in units_map:
                if trace.stats.location == ch_map['loc'] and trace.stats.channel == ch_map['channel']:
                    calc_output = output_dictionary[ch_map['unit'].lower()]
                    unit = ch_map['unit']
                    break
        else:
            calc_output = export_mode
            unit = unit_dictionary[calc_output]

        if conversion_map:
            for conv in conversion_map:
                if calc_output == conv["source"]:
                    conversion_factor = conv["factor"]
                    unit = conv["label"]
        else:
            conversion_factor = 1.0
            unit = unit_dictionary[calc_output]

        f, Pxx_den = fft_welch(trace, conversion_factor)

        new_psd = {
            "trace": trace,
            "f": f,
            "Pxx_den": Pxx_den,
            "unit": unit
        }
        psds.append(new_psd)

    return psds


def detrended_rms(traces, inv):
    import scipy.signal as sig
    import numpy as np

    for trace in traces:
        target = trace.copy()
        detrended = sig.detrend(target.data)
        rms = np.sqrt(np.mean(detrended**2))
        trace_name = target.stats.network + "." + target.stats.station + "." + target.stats.location + "." + target.stats.channel
        print(trace_name)
        print("RMS[Counts]: %f" % rms)

        rms_volts = (rms / 8388607) * 10
        print("RMS[Volts]: %f" % rms_volts)

        for i, network in enumerate(inv):
            if network.code == target.stats.network:
                net = i
                for j, station in enumerate(inv[net]):
                    if station.code == target.stats.station:
                        sta = j
                        for k, channel in enumerate(inv[net][sta]):
                            if channel.code == target.stats.channel and channel.location_code == target.stats.location:
                                sensitivity = channel.response.instrument_sensitivity.value

                                rms_ms2 = rms / sensitivity
                                print("RMS[m/s^2]: %f\n" % rms_ms2)
