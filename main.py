#!/usr/bin/env python3
import os
import sys
import time
import json

from datetime import datetime
from datetime import timedelta
# :units = "minutes since 1900-01-01 00:00:00";

date_time_str = '1900-01-01 00:00:00'
date_time_origin = datetime.strptime(date_time_str, '%Y-%m-%d %H:%M:%S')


# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

def load(path):

    with open(path, 'r') as fi:
        #print(fi.readlines())
        ply = fi.readlines()

    groups_ct = 0

    for line in ply:
        if line == '\n':
            groups_ct += 1

    mset = ply[0].split(',')
    print(f'path, {path} records: {len(mset)} groups: {len(ply)} {groups_ct}')

    return ply
    #     json_plist = json.load(fi)
    # self.VARS = json_plist


def build_from_time(time_strip_array):
    values = time_strip_array[0].split(',')
    time_header = []
    keys_header = {'records': 0}

    for t in values:
        xt = float(t.rstrip())
        delta_readable = date_time_origin + timedelta(minutes=xt)
        d = {'t': "%.2f" % xt, 'ts': str(delta_readable), 'data': []}
        time_header.append(d)

    keys_header['records'] = len(time_header)
    keys = ['lon', 'lat', 'depth']
    for i in keys:
        datum = load(i + '.csv')
        keys_header[i] = datum[0].rstrip().split(',')

    # depths = load('depth.csv')
    # depth_values = depths[0].split(',')
    group_len = len(keys_header['depth'])

    group_count = 0
    index = 0
    data = load('thetao.csv')

    grouped = [[[]] for i in range(len(time_header)+1)]

    data_max = len(data)

    for n, line in enumerate(data):
        if line == '\n' or n == data_max-1:
            group_count += 1
            if group_count == group_len:
                time_header[index]['data'] = grouped[index].copy()
                index += 1
                group_count = 0
            else:
                grouped[index].append([])
        else:
            grouped[index][group_count].append(line.rstrip().split(','))

    for k, header in enumerate(time_header):
        filename = f'data/data-{k+1}.json'
        with open(filename, 'w') as fp:
            json.dump(header, fp, indent=2)

        print(filename, '%ikb' % (os.path.getsize(filename)/1000), '')

    filename = f'data/data-keys.json'
    with open(filename, 'w') as fp:
        json.dump(keys_header, fp)  #;//, indent=2)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #//TODO: Dimensions: lon,lat,depth,time
    group = ['lon', 'lat', 'depth', 'time', 'thetao']

    build_from_time(load('time.csv'))

    #//TODO: build table keys for 'lon', 'lat', 'depth'


    # for i in group:
    #
    #     load(i+'.csv')
