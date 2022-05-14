#!/usr/bin/env python3

from json_encoder import JsonSafeEncoder
import marine_regions
import json

if __name__ == '__main__':

    elements = marine_regions.load()

    print(json.dumps(elements, indent=1, cls=JsonSafeEncoder))
