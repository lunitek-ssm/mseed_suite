from obspy import Stream
from obspy import Trace
from obspy.core.trace import Stats
from numpy import zeros


def reorder_stream(stream, channels_list):
	reordered_stream = Stream()
	for channel in channels_list:
		for trace in stream:
			tr_name = trace.stats.network + "." + trace.stats.station + "." + trace.stats.location + "." + trace.stats.channel
			if tr_name == channel:
				reordered_stream += trace

	# for i, channel in enumerate(channels_list):

		# print("STREAM OLD: \t" + stream[i].stats.network + "." + stream[i].stats.station + "." + stream[i].stats.location + "." + stream[i].stats.channel)
		# print("LIST: \t\t" + channel)
		# print("NEW STREAM: \t" + reordered_stream[i].stats.network + "." + reordered_stream[i].stats.station + "." + reordered_stream[i].stats.location + "." + reordered_stream[i].stats.channel)
		# print("")

	return reordered_stream


def find_trace(stream, channel):
	for i, trace in enumerate(stream):
		tr_name = trace.stats.network + "." + trace.stats.station + "." + trace.stats.location + "." + trace.stats.channel
		if tr_name == channel:
			return i
	return None


def fill_stream(stream, channels_list):
	for channel in channels_list:
		if find_trace(stream, channel) is None:
			data = zeros(len(stream[0]))
			stats = Stats()
			ch_parts = channel.split(".")
			stats.network = ch_parts[0]
			stats.station = ch_parts[1]
			stats.location = ch_parts[2]
			stats.channel = ch_parts[3]
			stats.starttime = stream[0].stats.starttime
			stats.sampling_rate = stream[0].stats.sampling_rate
			stats.delta = stream[0].stats.delta
			stats.npts = stream[0].stats.npts
			stats.calib = stream[0].stats.calib
			tr = Trace(data, stats)
			stream += tr
