import os
import glob
import zipfile
import math
from time import sleep

import pandas as pd
from pygbif import species
from pygbif import occurrences as occ
from pygbif.occurrences.download import GbifDownload

SLEEP_DURATION = 20

def getskey(z):
    return species.name_backbone(z)['usage']['key']

def run_gbif(my_vars, startyear, outfolder, today, userdata, skip=[]):

    gbif_obs = {}
    gbif_observers = {}

    for i, (k, v) in enumerate(my_vars.items()):

        if k in skip:
            continue

        records = []
        print("downloading ", k)

        splist = v['sci_name'].sum()
        spkeys = [getskey(x) for x in splist]
        spkeys = list(map(str, spkeys))

        gbif_query = GbifDownload(
            userdata.get('GBIF_USER'),
            userdata.get('GBIF_EMAIL')
        )

        gbif_query.add_predicate_dict({
            "type": "in",
            "key": "BASIS_OF_RECORD",
            "values": ['HUMAN_OBSERVATION', 'OBSERVATION', 'MACHINE_OBSERVATION', 'LIVING_SPECIMEN', 'MATERIAL_SAMPLE'],
            "matchCase": "false"
        })

        gbif_query.add_predicate_dict({
            "type": "equals",
            "key": 'HAS_COORDINATE',
            'value': 'TRUE',
            "matchCase": "false"
        })

        gbif_query.add_predicate_dict({
            "type": "equals",
            "key": 'HAS_GEOSPATIAL_ISSUE',
            'value': 'FALSE',
            "matchCase": "false"
        })

        gbif_query.add_predicate_dict({
            "type": "within",
            "geometry": "POLYGON((-100.551 36.917,-71.79 36.917,-71.79 49.612,-100.551 49.612,-100.551 36.917))"
        })

        gbif_query.add_predicate_dict({
            "type": "greaterThanOrEquals",
            "key": 'YEAR',
            'value': startyear,
            "matchCase": "false"
        })

        gbif_query.add_predicate_dict({
            "type": "in",
            "key": 'TAXON_KEY',
            "values": spkeys,
            "matchCase": "false"
        })

        xx = gbif_query.post_download(
            userdata.get('GBIF_USER'),
            userdata.get('GBIF_PWD')
        )

        while True:
            print(f"waiting to get download {xx}...")
            status = occ.download_meta(key=xx)['status']

            if status not in ['PREPARING', 'RUNNING']:
                if status == 'SUCCEEDED':
                    print("Download is ready, getting it")

                    output_path = k + "_gbif_obs"

                    if os.path.exists(output_path):
                        files = glob.glob(output_path + '/*.zip')
                        for f in files:
                            os.remove(f)
                    else:
                        os.mkdir(output_path)

                    occ.download_get(xx, output_path)

                else:
                    print("Status is", status)
                    print(occ.download_meta(key=xx))

                break

            sleep(SLEEP_DURATION)

        # --------------------------
        # read zip (unchanged logic)
        # --------------------------

        output_path = k + "_gbif_obs"
        files = os.listdir(output_path)
        file_path = os.path.join(output_path, files[0])

        print(file_path)

        base_name, extension = os.path.splitext(files[0])
        zf = zipfile.ZipFile(file_path)
        df = pd.read_csv(zf.open(base_name + '.csv'), sep='\t')

        df = df.loc[df['occurrenceStatus'] == "PRESENT"]

        print(len(df.index), k, " records")

        gbif_observers[k] = df['recordedBy'].unique()

        df = df[['gbifID', 'species', 'eventDate', 'decimalLatitude', 'decimalLongitude']]

        df.rename(columns={
            'gbifID': 'uid',
            "species": "sci_name",
            "eventDate": "obs_date",
            "decimalLatitude": "lat_dec",
            "decimalLongitude": "lon_dec"
        }, inplace=True)

        df["source"] = "GBIF"

        df['obs_date'] = df['obs_date'].str.slice(0, 10)
        df['obs_date'] = pd.to_datetime(df['obs_date'], format='mixed')

        df = df.dropna(subset='obs_date')
        df.reset_index(drop=True, inplace=True)

        df.to_csv(
            output_path + '/' + k + '_obs_gbifraw_' + today + '.csv',
            index=False
        )

        df_obs = pd.DataFrame(gbif_observers[k])
        df_obs.to_csv(
            output_path + '/' + k + '_observers_gbif_' + today + '.csv',
            index=False
        )

        gbif_obs[k] = df

        print("finished with ", k)

    return gbif_obs, gbif_observers


def run_gbif_pseudoabsences(my_vars, gbif_observers, startyear, userdata, skip=[]):

    from pygbif import species as species
    from pygbif import occurrences as occ
    from pygbif.occurrences.download import GbifDownload

    import os, math, zipfile
    from time import sleep
    import pandas as pd

    pabs = {}

    for i, (k, v) in enumerate(my_vars.items()):
        if k in skip:
            continue

        print("downloading ", k)

        observers = gbif_observers[k]
        observers = [x for x in observers if x == x]
        observers.sort(key=str.lower)

        start = 0
        end = len(observers)
        step = 50

        for x in range(start, end, step):
            observers_chunk = observers[x:x+step]

            print(f'chunk {x}: ', observers_chunk[0])

            gbif_query = GbifDownload(
                userdata.get('GBIF_USER'),
                userdata.get('GBIF_EMAIL')
            )

            gbif_query.add_predicate_dict({"type": "equals", "key": 'HAS_COORDINATE', 'value': 'TRUE', "matchCase": "false"})
            gbif_query.add_predicate_dict({"type": "equals", "key": 'HAS_GEOSPATIAL_ISSUE', 'value': 'FALSE', "matchCase": "false"})
            gbif_query.add_predicate_dict({"type": "within", "geometry": "POLYGON((-100.551 36.917,-71.79 36.917,-71.79 49.612,-100.551 49.612,-100.551 36.917))"})
            gbif_query.add_predicate_dict({"type": "greaterThanOrEquals", "key": 'YEAR', 'value': startyear, "matchCase": "false"})
            gbif_query.add_predicate_dict({"type": "in", "key": 'RECORDED_BY', 'values': observers_chunk, "matchCase": "false"})

            # taxon filters (unchanged)
            if k == 'fish':
                gbif_query.add_predicate_dict({"type": "equals", "key": 'TAXON_KEY', 'value': 44, "matchCase": "false"})
            if k == 'plant':
                gbif_query.add_predicate_dict({"type": "equals", "key": 'TAXON_KEY', 'value': 6, "matchCase": "false"})
            if k == 'invert':
                gbif_query.add_predicate_dict({"type": "equals", "key": 'TAXON_KEY', 'value': 1, "matchCase": "false"})
                gbif_query.add_predicate_dict({"type": "not", "predicate": {"type": "equals", "key": 'TAXON_KEY', 'value': 44}})

            xx = gbif_query.post_download(
                userdata.get('GBIF_USER'),
                userdata.get('GBIF_PWD')
            )

            while True:
                print(f"waiting to get download {xx}...")
                status = occ.download_meta(key=xx)['status']

                if status not in ['PREPARING', 'RUNNING']:
                    if status == 'SUCCEEDED':
                        output_path = k+"_gbif_observer_potential_pseudoabs"
                        if not os.path.exists(output_path):
                            os.mkdir(output_path)

                        occ.download_get(xx, output_path)
                    break

                sleep(20)

    return pabs
