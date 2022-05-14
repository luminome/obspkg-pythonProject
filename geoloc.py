# -*- coding: utf-8 -*-

from geopy.geocoders import Nominatim
import json
import asyncio

# create a Nominatim object
geolocator = Nominatim(user_agent="sac_app_test")

sanitized = [
    'place_id',
    'boundingbox',
    'display_name',
    'lon',
    'lat',
    'type',
    'extratags'
]

locations_group = [
    'marseille',
    'rome',
    'toulon',
    'tarragona',
    'palermo',
    'beirut',
    'nice',
    'cannes',
    'genoa',
    'livorno',
    'athens',
    'antalia',
    'alexandria',
    'gibraltar',
    'ajaccio',
    'algiers',
    'tripoli',
    'venice',
    'trieste',
    'zadar',
    'split',
    'dubrovnik',
    'durres'
]

locations_group.sort()
print(locations_group)

loc_index = 0
locations_all = {}

loop = asyncio.get_event_loop()


def recover(name):
    # Opening JSON file
    filename = f'data/data-{name}-test.json'
    f = open(filename)
    data = json.load(f)
    return data


def deliver(data, name):
    filename = f'data/data-{name}-test.json'
    with open(filename, 'w') as fp:
        json.dump(data, fp, indent=2)


def shutdown():
    for task in asyncio.all_tasks():
        task.cancel()


def callback():
    global loc_index

    while True:
        try:
            place = locations_group[loc_index]
            loc_index += 1
            if place.title() not in locations_all:
                print(place, "loading request")
                break
            else:
                print(place, "skipped")

        except IndexError:
            shutdown()
            print('end of list reached')
            return

    result = {}

    location = geolocator.geocode(place, language='en', namedetails=True, addressdetails=True, extratags=True, featuretype='city')
    for n in sanitized:
        result[n] = location.raw[n] if n != 'display_name' else location.raw[n]  #;//.encode('utf-8').decode('utf-8')

    result['country'] = location.raw['address']['country']

    locations_all[place.title()] = result

    if loc_index < len(locations_group):
        loop.call_later(2, callback)
    else:
        shutdown()


async def main():
    while True:
        await asyncio.sleep(1)


if __name__ == '__main__':

    locations_all = recover('locs')

    try:
        loop.call_later(1, callback)
        loop.run_until_complete(main())
    except asyncio.exceptions.CancelledError:

        deliver(locations_all, 'locs')
        exit()

