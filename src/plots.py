from math import pi
import numpy
import matplotlib.pyplot as plt
from obspy.imaging.cm import pqlx
from obspy.signal.spectral_estimation import get_nlnm as get_NHNM
from obspy.signal.spectral_estimation import get_nhnm as get_NLNM

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
    "INFRA": "PA"
}

color_dictionary = {
    "ACC": "red",
    "VEL": "blue",
    "DISP": "green",
    "INFRA": "purple",
    "NATURAL": "orange"
}


color_infrasound_array = ["red", "blue", "green", "purple", "orange", "grey"]

def raw(traces, start_t, path, print_output=False):
    filename = path + start_t.strftime("%Y%m%d%H%M%S") + "_traces_raw.pdf"
    fig = traces.plot(handle=True)
    for ax in fig.axes:
        ax.set_ylabel("Counts")
    plt.savefig(filename, bbox_inches='tight')
    plt.cla()
    plt.clf()
    plt.close()
    if print_output:
        print("Plotted raw traces")


def dayplot(traces, start_t, path, print_output=False):
    for tr in traces:
        tr.plot(outfile=path + start_t.strftime("%Y%m%d%H%M%S") + "-" +
        tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + ".dayplot.pdf", type="dayplot")
    if print_output:
        print("Plotted traces dayplot")


def deconvoluted(traces, units_map, mode, start_t, path, print_output=False):

    for trace in traces:
        if mode == "NATURAL":
            for ch_map in units_map:
                if trace.stats.location == ch_map['loc'] and trace.stats.channel == ch_map['channel']:
                    calc_output = output_dictionary[ch_map['unit'].lower()]
                    unit = ch_map['unit']
                    break
        else:
            calc_output = mode
            unit = unit_dictionary[calc_output]

        outfile = path + start_t.strftime("%Y%m%d%H%M%S")
        outfile += "_" + trace.stats.network + "." + trace.stats.station + "." + trace.stats.location + "." + trace.stats.channel + "_"
        outfile += calc_output + ".pdf"

        color = color_dictionary[calc_output]
        label = unit
        # label = "Volts"
        fig = trace.plot(handle=True, color=color)
        for ax in fig.axes:
            ax.set_ylabel(label)
            ax.yaxis.set_label_text(label)
        plt.savefig(outfile, bbox_inches='tight')
        plt.cla()
        plt.clf()
        plt.close()

    if print_output:
        print("Plotted deconvoluted traces")


def ppsd(ppsds, structural_mode, path, print_output=False):
    ERROR = False
    for ppsd in ppsds:
        ppsd["data"].save_npz(path + "npz/" + ppsd["name"] + ".npz")
        try:
            if structural_mode:
                ppsd["data"].plot(cmap=pqlx, filename=path + "plot/" + ppsd["name"] + ".pdf", period_lim=(0.02, 10), show_mean=True, show_noise_models=False)
            else:
                ppsd["data"].plot(cmap=pqlx, filename=path + "plot/" + ppsd["name"] + ".pdf", show_mean=True, show_noise_models=True)

        except Exception as e:
            print("ERROR: Cannot plot PPSD - %s" % str(e))
            ERROR = True
    return ERROR


def psd(psds, units_map, mode, start_t, path):
    for psd in psds:

        psd_file_name = start_t.strftime("%Y%m%d%H%M%S") + '.'
        psd_file_name += psd["trace"].stats.network + '.' + psd["trace"].stats.station + '.' + psd["trace"].stats.location + "." + psd["trace"].stats.channel + ".psd"
        outfile = path + psd_file_name + ".pdf"

        if mode == "NATURAL":
            for ch_map in units_map:
                if psd["trace"].stats.location == ch_map['loc'] and psd["trace"].stats.channel == ch_map['channel']:
                    calc_output = output_dictionary[ch_map['unit'].lower()]
                    break
        else:
            calc_output = mode

        color = color_dictionary[calc_output]

        f_zero = psd["f"][1:]
        Pxx_den_zero = psd["Pxx_den"][1:]

        plt.autoscale(axis='y')
        plt.semilogx(f_zero, Pxx_den_zero, label=psd["trace"].stats.location + " - " + psd["trace"].stats.channel, color=color)
        plt.legend()

        plt.xlabel('Frequency [Hz]')
        plt.ylabel('PSD [(' + psd["unit"] + ')^2/Hz]')

        plt.savefig(outfile, bbox_inches='tight')
        plt.clf()
        plt.cla()
        plt.close()


def psd_seismic(psds, units_map, mode, start_t, path):
    for psd in psds:

        psd_file_name = start_t.strftime("%Y%m%d%H%M%S") + '.'
        psd_file_name += psd["trace"].stats.network + '.' + psd["trace"].stats.station + '.' + psd["trace"].stats.location + "." + psd["trace"].stats.channel + ".psd"
        outfile = path + psd_file_name + ".pdf"

        if mode == "NATURAL":
            for ch_map in units_map:
                if psd["trace"].stats.location == ch_map['loc'] and psd["trace"].stats.channel == ch_map['channel']:
                    calc_output = output_dictionary[ch_map['unit'].lower()]
                    break
        else:
            calc_output = mode

        color = color_dictionary[calc_output]

        f_zero = psd["f"][1:]
        Pxx_den_zero = psd["Pxx_den"][1:]

        plt.autoscale(axis='y')

        NLNMper, NLNMpower = get_NLNM()
        NHNMper, NHNMpower = get_NHNM()

        if psd["unit"] == "ACC":
            plt.plot(NLNMper, NLNMpower, 'k')
            plt.plot(NHNMper, NHNMpower, 'k')

        # Velocità o spostamento
        if psd["unit"] == "VEL" or psd["unit"] == "DISP":
            nl_start_vector = 2 * numpy.log10(NLNMper)
            nl_conv_factor = numpy.log10(2) + numpy.log10(pi)
            nl_scalar = numpy.full(nl_start_vector.shape, nl_conv_factor)
            nl_conv_vector = nl_start_vector + nl_scalar
            NLNMpower += nl_conv_vector

            nh_start_vector = 2 * numpy.log10(NHNMper)
            nh_conv_factor = numpy.log10(2) + numpy.log10(pi)
            nh_scalar = numpy.full(nh_start_vector.shape, nh_conv_factor)
            nh_conv_vector = nh_start_vector + nh_scalar

            NHNMpower += nh_conv_vector

            if psd["unit"] == "VEL":
                plt.plot(NLNMper, NLNMpower, 'k')
                plt.plot(NHNMper, NHNMpower, 'k')

            # psd_ascii(NLNMper, NLNMpower, "conv1_NLNMper.txt")
            # psd_ascii(NHNMper, NHNMpower, "conv1_NHNMper.txt")

        # Spostamento
        if psd["unit"] == "DISP":
            nl_start_vector = 2 * numpy.log10(NLNMper)
            nl_conv_factor = numpy.log10(2) + numpy.log10(pi)
            nl_scalar = numpy.full(nl_start_vector.shape, nl_conv_factor)
            nl_conv_vector = nl_start_vector + nl_scalar
            NLNMpower += nl_conv_vector

            nh_start_vector = 2 * numpy.log10(NHNMper)
            nh_conv_factor = numpy.log10(2) + numpy.log10(pi)
            nh_scalar = numpy.full(nh_start_vector.shape, nh_conv_factor)
            nh_conv_vector = nh_start_vector + nh_scalar

            NHNMpower += nh_conv_vector

            plt.plot(NLNMper, NLNMpower, 'k')
            plt.plot(NHNMper, NHNMpower, 'k')

        # Infrasuono

        if psd["unit"] == "PA":
            GLNper = [0.143, 0.167, 0.2, 0.25, 0.333, 0.5, 1.0, 1.111, 1.25, 1.429, 1.667, 2.0, 2.5, 3.333, 5.0, 6.667, 10.0, 11.765 ,14.286, 20.0, 25.0, 33.333]
            GLNpower = [-82.0, -80.0, -75.0, -72.0, -71.0, -68.0, -64.0, -62.0, -60.0, -59.0, -57.0, -55.0, -51.0, -47.0, -40.0, -42.0, -47.0, -44.0, -41.0, -36.0, -33.0, -29.0]
            GLWper = [0.143, 0.167, 0.2, 0.25, 0.333, 0.5, 1.0, 1.111, 1.25, 1.429, 1.667, 2.0, 2.5, 3.333, 5.0, 6.667, 10.0, 11.765 ,14.286, 20.0, 25.0, 33.333]
            GLWpower = [-61.0, -61.0, -61.0, -60.0, -58.0, -56.0, -53.0, -51.0, -50.0, -48.0, -45.0, -43.0, -40.0, -35.0, -28.0, -31.0, -33.0, -31.0, -30.0, -25.0, -20.0, -18.0]
            plt.plot(GLNper, GLNpower, 'k')
            plt.plot(GLWper, GLWpower, 'k')            


        plt.semilogx(1 / f_zero, 10 * numpy.log10(abs(Pxx_den_zero)), label="PSD " + psd["trace"].stats.location + " - " + psd["trace"].stats.channel, color=color)
        # CONFRONTO CON TEST DI SLEEMAN PER DEBUG
        # plt.plot(1 / fre, 10 * numpy.log10(abs(cpval)), 'r', label="Sleeman " + trace.stats.station + ' ' + trace.stats.location + ' ' + trace.stats.channel)
        plt.legend()
        plt.grid(which="both")
        plt.xlim((numpy.amin(1 / f_zero), numpy.amax(1 / f_zero)))
        # plt.xlim((numpy.amin(1 / fre), numpy.amax(1 / fre)))
        # plt.xlim((0.5, 20))
        plt.xlabel('Period [S]')
        plt.ylabel('Seismic PSD [(' + psd["unit"] + ')^2/Hz]')
        plt.savefig(outfile, bbox_inches='tight')
        plt.clf()
        plt.cla()
        plt.close()


def psd_infrasound(psds, units_map, mode, start_t, path):

    psd = psds[0]

    psd_file_name = start_t.strftime("%Y%m%d%H%M%S") + '.'
    psd_file_name += psd["trace"].stats.network + '.' + psd["trace"].stats.station + ".psd"
    outfile = path + psd_file_name + ".pdf"

    psd_index = 0

    plt.ylim(-100.0, 0.0)

    if psd["unit"] == "PA":
        GLNper = [0.143, 0.167, 0.2, 0.25, 0.333, 0.5, 1.0, 1.111, 1.25, 1.429, 1.667, 2.0, 2.5, 3.333, 5.0, 6.667, 10.0, 11.765 ,14.286, 20.0, 25.0, 33.333]
        GLNpower = [-82.0, -80.0, -75.0, -72.0, -71.0, -68.0, -64.0, -62.0, -60.0, -59.0, -57.0, -55.0, -51.0, -47.0, -40.0, -42.0, -47.0, -44.0, -41.0, -36.0, -33.0, -29.0]
        GLWper = [0.143, 0.167, 0.2, 0.25, 0.333, 0.5, 1.0, 1.111, 1.25, 1.429, 1.667, 2.0, 2.5, 3.333, 5.0, 6.667, 10.0, 11.765 ,14.286, 20.0, 25.0, 33.333]
        GLWpower = [-61.0, -61.0, -61.0, -60.0, -58.0, -56.0, -53.0, -51.0, -50.0, -48.0, -45.0, -43.0, -40.0, -35.0, -28.0, -31.0, -33.0, -31.0, -30.0, -25.0, -20.0, -18.0]
        plt.plot(GLNper, GLNpower, 'k')
        plt.plot(GLWper, GLWpower, 'k')      



    for psd in psds:

        color = color_infrasound_array[psd_index%6]
        psd_index += 1

        f_zero = psd["f"][1:]
        Pxx_den_zero = psd["Pxx_den"][1:]




        # Infrasuono

      


        plt.semilogx(1 / f_zero, 10 * numpy.log10(abs(Pxx_den_zero)), label="PSD " + psd["trace"].stats.location + " - " + psd["trace"].stats.channel, color=color)
        # CONFRONTO CON TEST DI SLEEMAN PER DEBUG
        # plt.plot(1 / fre, 10 * numpy.log10(abs(cpval)), 'r', label="Sleeman " + trace.stats.station + ' ' + trace.stats.location + ' ' + trace.stats.channel)
    plt.legend()
    plt.grid(which="both")
    plt.xlim((numpy.amin(1 / f_zero), numpy.amax(1 / f_zero)))
        # plt.xlim((numpy.amin(1 / fre), numpy.amax(1 / fre)))
        # plt.xlim((0.5, 20))
    plt.xlabel('Period [S]')
    plt.ylabel('Seismic PSD [(' + psd["unit"] + ')^2/Hz]')
    plt.savefig(outfile, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()



def spectrogram(traces, path, print_output):
    for tr in traces:
        tr.spectrogram(outfile=path + tr.stats.starttime.strftime("%Y%m%d%H%M%S") + "-" +
        tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + ".spectrogram.pdf")
    if print_output:
        print("Calculated and plotted spectrogram")
