import os
import xml.dom.minidom
import obspy
from obspy.clients.nrl import NRL
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.core.inventory.util import Equipment


def explore(version, keys_list, nrl):
    if keys_list[-1] == "FOUND":
        return keys_list
    if keys_list[-1] == "NOT_FOUND":
        return ["{\"error\":true}"]

    keys_str = ""
    for key in keys_list:
        keys_str += "[\"" + key + "\"]"

    # print(keys_str)
    found = eval("nrl.sensors" + keys_str)
    if type(found) is tuple:
        if found[0] == version:
            keys_list.append("FOUND")
        elif 'RESP' in found[1]:
            keys_list.pop()
            keys_list.append("NOT_FOUND")
        else:
            keys_list.append(found[0])
    else:
        for key, val in found.items():
            keys_list.append(key)
            keys_list = explore(version, keys_list, nrl)
            if keys_list[-1] == "FOUND":
                break
            else:
                keys_list.pop()

    return keys_list


def get_nrl_catalog(iris_path):
    return NRL(os.path.abspath(iris_path))


def get_element_index(inv_list, name):

    for i, elem in enumerate(inv_list):
        if elem["code"] == name:
            return i
    return None


def stationxml_from_metadata(pc_file, fd_file, nrl):
    chs_pc = []
    chs_fd = []
    channels_map = []

    with open(fd_file) as factory_data_xml:
        DOMTree_fd = xml.dom.minidom.parse(factory_data_xml)
        doc_fd = DOMTree_fd.documentElement
        device_fd = doc_fd.getElementsByTagName("device")
        model_fd = device_fd[0].getAttribute("model")
        ss4_fd = doc_fd.getElementsByTagName("ss4")
        channels_fd = ss4_fd[0].getElementsByTagName("channel")

        for channel in channels_fd:
            ch = {}
            ch['index'] = channel.getAttribute("index")
            ch['model'] = model_fd
            if channel.hasAttribute("hwfs"):
                ch['hwfs'] = channel.getAttribute("hwfs")
            else:
                ch['hwfs'] = ss4_fd[0].getAttribute("hwfs")
            chs_fd.append(ch)

    with open(pc_file) as poseidon_config_xml:
        DOMTree_pc = xml.dom.minidom.parse(poseidon_config_xml)
        doc_pc = DOMTree_pc.documentElement
        boards = doc_pc.getElementsByTagName("board")
        if not boards:
            boards = doc_pc.getElementsByTagName("sensor")
        samplerate_tags = doc_pc.getElementsByTagName("samplerate")
        samplerate = samplerate_tags[0].getAttribute("freq1")

        # READ PRETRIGGER VALUE
        trigger_tag = doc_pc.getElementsByTagName("trigger")
        pretrigger = int(trigger_tag[0].getAttribute("pretrigger"))

        for board in boards:
            channels = board.getElementsByTagName("channel")
            for channel in channels:
                ch = {}
                miniseed_ch = channel.getElementsByTagName("miniseed_ch")

                ch['net'] = miniseed_ch[0].getAttribute("net")
                ch['sta'] = miniseed_ch[0].getAttribute("sta")
                ch['name'] = miniseed_ch[0].getAttribute("cha")
                ch['location'] = miniseed_ch[0].getAttribute("loc")
                ch['pga'] = channel.getAttribute("pga")
                ch['filt'] = channel.getAttribute("filt")
                ch['index'] = channel.getAttribute("index")
                ch['sample_rate'] = samplerate

                channel_convert = {
                    "loc": ch["location"],
                    "channel": ch["name"],
                    "name": board.getAttribute("name")
                }
                channels_map.append(channel_convert)

                # GET SENSOR SEARCH KEYS
                producer = board.getAttribute("producer")
                model = board.getAttribute("model")
                version = board.getAttribute("version")
                sentype = board.getAttribute("type")
                ch['producer'] = producer
                ch['version'] = version
                ch['type'] = sentype
                keys_list = [producer, model]
                keys = explore(version, keys_list, nrl)
                keys.pop()
                ch['keys'] = keys

                chs_pc.append(ch)

    # print(CHANNELS_MAP)

    inventory = {
        "networks": [],
        "source": "Lunitek Mseed Suite"
    }

    network = {
        "code": chs_pc[0]['net'],
        "description": "",
        "stations": []
    }

    station = {
        "code": chs_pc[0]['sta'],
        "elevation": 0.0,
        "site": "",
        "longitude": 0.0,
        "latitude": 0.0,
        "channels": []
    }

    for i in range(len(chs_pc)):

        net_idx = get_element_index(inventory['networks'], chs_pc[i]['net'])
        if net_idx is None:
            tmp_net = {
                "code": chs_pc[i]['net'],
                "description": "",
                "stations": []
            }
            inventory['networks'].append(tmp_net)
            net_idx = len(inventory['networks']) - 1

        sta_idx = get_element_index(inventory['networks'][net_idx]['stations'], chs_pc[i]['sta'])
        if sta_idx is None:
            tmp_sta = {
                "code": chs_pc[i]['sta'],
                "elevation": 0.0,
                "site": "",
                "longitude": 0.0,
                "latitude": 0.0,
                "channels": []
            }
            inventory['networks'][net_idx]['stations'].append(tmp_sta)
            sta_idx = len(inventory['networks'][net_idx]['stations']) - 1

        if "Atlas" in chs_fd[i]['model'] or "Triton" in chs_fd[i]['model'] or "MCDR" in chs_fd[i]['model']:
            vpp = (int(chs_fd[i]['hwfs']) * 2) / int(chs_pc[i]['pga'])
            vpp_str = str(int(vpp)) + " Vpp"
            filter_name = "Minimum"
            if chs_pc[i]['filt'] == "LIN":
                filter_name = "Linear"
            datalogger_keys = ['Lunitek', 'ATLAS', vpp_str, samplerate, filter_name]
            datalogger_name = "ATLAS, gain 1, " + vpp_str + ", " + samplerate + " sps, " + filter_name + " Phase filter"

        elif "Sentinel" in chs_fd[i]["model"]:
            if int(chs_fd[i]['index']) < 3:
                datalogger_keys = ['Lunitek', 'Sentinel', "MEMS Accelerometer (10 Vpp)", samplerate]
                datalogger_name = "Sentinel r.1, gain 1, 10 Vpp, " + samplerate + " sps"
            else:
                vpp = chs_fd[i]['hwfs']
                vpp_str = vpp + " Vpp"
                datalogger_keys = ['Lunitek', 'Sentinel', "Geophone", vpp_str, samplerate]
                datalogger_name = "Sentinel r.1, gain 1, " + vpp_str + " , " + samplerate + " sps"

        tmp_cha = {
            "code": chs_pc[i]['name'],
            "elevation": 0.0,
            "depth": 0.0,
            "longitude": 0.0,
            "location_code": chs_pc[i]['location'],
            "sample_rate": samplerate,
            "azimuth": 0.0,
            "latitude": 0.0,
            "dip": 0.0,
            "sensor": chs_pc[i]['keys'],
            "datalogger": datalogger_keys,
            "sensor_name": chs_pc[i]['version'],
            "sensor_type": chs_pc[i]['type'],
            "sensor_manufacturer": chs_pc[i]['producer'],
            "datalogger_name": datalogger_name
        }

        inventory['networks'][net_idx]['stations'][sta_idx]['channels'].append(tmp_cha)
        # station['channels'].append(tmp_cha)

    # network['stations'].append(station)
    # inventory['networks'].append(network)

    # print(inventory)

    stationxmlfile = network['code'] + "." + station['code'] + ".station.xml"

    try:
        inv = Inventory(
            networks=[],
            source=inventory['source'])

        for net in inventory['networks']:
            new_network = Network(
                code=net['code'],
                description=net['description'],
                start_date="2000-01-01T00:00:00",
                stations=[]
            )
            for sta in net['stations']:
                new_station = Station(
                    code=sta['code'],
                    latitude=sta['latitude'],
                    longitude=sta['longitude'],
                    elevation=sta['elevation'],
                    creation_date=obspy.UTCDateTime.now(),
                    site=Site(name=sta['site']),
                    channels=[]
                )
                for cha in sta['channels']:
                    # print("setting sensor description")
                    sensor = Equipment(type=cha['sensor_type'],
                                    description=cha['sensor_name'],
                                    manufacturer=cha['sensor_manufacturer'],
                                    vendor=cha['sensor_manufacturer'],
                                    model=cha['sensor_name'],
                                    resource_id="Lunitek station XML builder")

                    datalogger = Equipment(type='datalogger',
                                    description=cha['datalogger_name'],
                                    manufacturer='Lunitek',
                                    vendor='Lunitek',
                                    model=cha['sensor_name'],
                                    resource_id="Lunitek station XML builder")
                    new_channel = Channel(
                        code=cha['code'],
                        location_code=cha['location_code'],
                        latitude=cha['latitude'],
                        longitude=cha['longitude'],
                        elevation=cha['elevation'],
                        depth=cha['depth'],
                        azimuth=cha['azimuth'],
                        dip=cha['dip'],
                        sample_rate=cha['sample_rate'],
                        sensor=sensor,
                        data_logger=datalogger
                    )

                    if cha['datalogger'][0] != "Lunitek":
                        cha['datalogger'].insert(0, "Lunitek")
                    # print(cha['sensor'])
                    # print(cha['datalogger'])
                    try:
                        response = nrl.get_response(
                            sensor_keys=cha['sensor'],
                            datalogger_keys=cha['datalogger']
                        )
                    except KeyError:
                        print("Could not find a sensors or datalogger with specified keys, exiting now")
                        print(KeyError)

                    new_channel.response = response

                    new_station.channels.append(new_channel)
                new_network.stations.append(new_station)
            inv.networks.append(new_network)
            inv.write(stationxmlfile, format="stationxml", validate=True)
        if not os.path.isfile(stationxmlfile):
            print("ERROR: Could not export StationXML File")
            exit(-1)
            # print("INFO: Succesfully created StationXML File")
        else:
            return stationxmlfile, channels_map, pretrigger
    except Exception as ex:
        print(ex)
        exit(-1)
