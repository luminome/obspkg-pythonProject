#!/usr/bin/env python3
import json
from typing import List
import numpy as np
import math

from shapely.geometry import Polygon, MultiPolygon, LinearRing, LineString, MultiLineString, box
from shapely import ops
from shapely.validation import make_valid

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection
from json_encoder import JsonSafeEncoder


def gpd_scope_shp_to_bounds(shp_path, bounds_box):
    import geopandas as gpd
    # minx, miny, maxx, maxy = bounds
    # geographic_bounds = box(minx, miny, maxx, maxy)

    regions = gpd.read_file(shp_path)
    geo_b = gpd.GeoDataFrame(index=[0], crs=regions.crs, geometry=[bounds_box])
    intersections2 = gpd.overlay(geo_b, regions, how='intersection', keep_geom_type=False)
    # intersections2.to_file("../data/goas_med/med_zone.shp")
    return intersections2


def show_progress(item_name, count, total):
    if total > 0:
        pit = math.ceil(((count) / (total-1)) * 100)
    else:
        pit = 'NaN'
    print(f"\r{item_name} {pit}% complete", end='')


def plot_polygon(ax, poly, **kwargs):
    path = Path.make_compound_path(
        Path(np.asarray(poly.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])

    patch = PathPatch(path, **kwargs)
    collection = PatchCollection([patch], **kwargs)

    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection


def package_poly(poly):
    poly_lex = {'out': coords_to_flat(poly.exterior.coords)}
    if len(poly.interiors):
        poly_lex['ins'] = [coords_to_flat(p.coords) for p in poly.interiors]
    return poly_lex


def refine_bounds(bounds, cfg, no_pad=False):
    tb = []
    print(bounds)

    for i, b in enumerate(bounds):
        if no_pad:
            tb.append(round(b))
        else:
            tb.append(math.floor(b) + cfg['bounds_keys'][i])

    w = abs(tb[0] - tb[2])
    h = abs(tb[1] - tb[3])
    return [tb, w, h]


def coords_to_flat(coords_list):
    coords_points = []
    list(coords_points.extend((round(cx, 4), round(cy, 4))) for cx, cy in coords_list)
    return coords_points


def get_data_from_sector(sector_tuple, sector_width, data_packet, level):
    # data, lons_ax, lats_ax, resolution, el_w)

    interval = data_packet['degs']
    data = data_packet['data']
    resolution = data_packet['reso'][level]
    lons_ax = data_packet['lons']
    lats_ax = data_packet['lats']
    indices = data_packet['indx']

    mask_lon = (lons_ax >= sector_tuple[0]) & (lons_ax < sector_tuple[0] + sector_width)
    mask_lat = (lats_ax <= sector_tuple[1]) & (lats_ax > sector_tuple[1] - sector_width)


    # try:
    xe = [np.where(lons_ax >= i)[0][0] for n, i in enumerate(lons_ax[mask_lon]) if n % resolution == 0]
    ye = [np.where(lats_ax <= i)[0][0] for n, i in enumerate(lats_ax[mask_lat]) if n % resolution == 0]
    # print(xe, ye)
    # except IndexError:
    #     return None

    if len(xe) == 0 or len(ye) == 0:
        return None

    xm = (np.amin(lons_ax[mask_lon])-sector_tuple[0])
    ym = sector_tuple[1]-(np.amax(lats_ax[mask_lat]))

    xxs = np.sign(xm)
    xxi = int(math.floor(abs(xm) * (interval/resolution))*xxs)
    yyi = int(math.floor(ym * (interval/resolution)))

    try:
        data_min = np.nanmin(data)
        data_max = np.nanmax(data)

    except TypeError:
        data_min = np.nan
        data_max = np.nan
        pass

    # mt = np.mean(data)

    data_master = []
    for n_i in range(0, len(indices[0])):
        if len(indices[1]) > 1:
            data_master = [[]] * len(indices[0])

        for n_j in range(0, len(indices[1])):
            #print(data.dtype)

            if data.dtype != '<U17':
                sector_data = np.zeros([int((interval*sector_width)/resolution), int((interval*sector_width)/resolution)])
                sector_data[:] = np.nan
            else:
                sector_data = np.empty([int((interval*sector_width)/resolution), int((interval*sector_width)/resolution)], dtype='<U17')

            if len(indices[0]) > 1 and len(indices[1]) == 1:
                # assume time data here
                for iy, dy in enumerate(ye):
                    for ix, dx in enumerate(xe):
                        try:
                            sector_data[iy+yyi, ix+xxi] = data[n_i, dy, dx]
                        except IndexError:
                            pass

                # print(data.dtype)
                # print(sector_data)

                # #// WHAT THE HELL IS THIS?
                # ty = np.where(np.isnan(sector_data), None, sector_data)
                # data_master.append(ty.tolist() if data.dtype == '<U17' else ty)
                data_master.append(sector_data.tolist())  #;//.tolist())


            elif len(indices[0]) > 1 and len(indices[1]) > 1:
                if n_j % 4 == 0:
                    # assume time data here
                    for iy, dy in enumerate(ye):
                        for ix, dx in enumerate(xe):
                            try:
                                sector_data[iy+yyi, ix+xxi] = data[n_i, n_j, dy, dx]
                            except IndexError:
                                pass

                    if not np.all(np.isnan(sector_data)):
                        ty = np.where(np.isnan(sector_data), None, sector_data)
                        data_master[n_i].append({str(indices[1][n_j]): ty.tolist()})

            else:
                for iy, dy in enumerate(ye):
                    for ix, dx in enumerate(xe):
                        try:
                            sector_data[iy+yyi, ix+xxi] = data[dy, dx]
                        except IndexError:
                            pass
                data_master = sector_data

    extents = [sector_tuple[0], sector_tuple[0]+sector_width, sector_tuple[1]-sector_width, sector_tuple[1]]
    return {'data': data_master, 'extents': extents, 'limits': [data_min, data_max]}


def prepare_map_points(the_poly, simplification, min_area):
    shapes = []
    interiors = []
    exterior = []
    multipoly = []
    group = None

    def qualify(ct, sub_shape, min_a):
        area_test = Polygon(sub_shape.coords)
        shape_simple = sub_shape.simplify(simplification, True)
        if area_test.area >= min_a:
            return [ct, shape_simple]
        return None

    if the_poly.type == 'Polygon':
        group = [the_poly]
    if the_poly.type == 'MultiPolygon':
        group = the_poly.geoms

    for the_poly in group:
        sub_interiors = []
        for c, n in enumerate(the_poly.interiors):
            e = qualify(c, n, min_area)
            if e is not None:
                sub_interiors.append(e[1])
                interiors.append(e[1])
                e.append('interior')
                shapes.append(e)

        e = qualify(0, the_poly.exterior, min_area)
        if e is not None:
            e.append('exterior')
            shapes.append(e)
            exterior.append(e[1])

            if len(the_poly.interiors):
                multipoly.append(Polygon(e[1], sub_interiors))
            else:
                multipoly.append(Polygon(e[1]))

    if len(group) == 1:
        print('prepare_map_points poly')
        return shapes, Polygon(exterior[0], interiors)
    else:
        print('prepare_map_points multipoly')
        return shapes, MultiPolygon(multipoly)


def get_map_shape_guides(master_poly, test_bounds, simps, areas, plt, buffer_constant=0.1, simplification_const=0.005):
    # #//TODO, make canonical shapes guide
    master_shapes, simplified_poly = prepare_map_points(master_poly, simps[3], areas[0])
    guides = []

    print("guide master_shapes", len(master_shapes))

    for c, n in enumerate(master_shapes):
        it = math.floor(((c) / len(master_shapes)) * 100)
        print("\rGuides {}% complete".format(it), end='')

        poly_c = test_bounds.intersection(Polygon(n[1]))
        poly_shape = test_bounds.difference(poly_c)

        if poly_shape.type == 'MultiPolygon':
            kbf = poly_shape.buffer(buffer_constant)
            kbf = kbf.simplify(simplification_const)
            print(kbf)

            if kbf.type == 'Polygon':
                poly_shape = kbf
            else:
                for w, pl in enumerate(kbf.geoms):
                    rtl = test_bounds.intersection(pl.exterior)
                    if rtl.type == 'MultiLineString':
                        dz = ops.linemerge(rtl)
                        print(dz.type)
                        if dz.type == 'LineString':
                            plt.plot(*dz.coords.xy, color="red", linewidth="4")
                            guides.append([n[0], len(guides), dz, None, dz.is_ring, n[1].length])
                        else:
                            for dzi in dz.geoms:
                                plt.plot(*dzi.coords.xy, color="red", linewidth="4")
                                guides.append([n[0], len(guides), dzi, None, dzi.is_ring, n[1].length])
                    else:
                        plt.plot(*rtl.coords.xy, color="blue", linewidth="4")
                        guides.append([n[0], len(guides), rtl, None, rtl.is_ring, n[1].length])

        if poly_shape.type == 'Polygon':
            if len(poly_shape.interiors):
                for i in poly_shape.interiors:
                    bf = buffer_constant
                    if n[2] == 'exterior':
                        bf *= -1

                    i_poly = Polygon(i)
                    kbf = i_poly.buffer(bf)
                    kbf = kbf.simplify(simplification_const)

                    if kbf.type == 'MultiPolygon' and n[2] == 'exterior':
                        s_poly = [sp for sp in kbf.geoms if sp.area > 1.0][0]
                        plt.plot(*s_poly.exterior.xy, color="green", linewidth="2")
                        guides.append([n[0], len(guides), s_poly.exterior, n[2],  True, n[1].length])
                    else:
                        plt.plot(*kbf.exterior.xy, color="black", linewidth="2")
                        guides.append([n[0], len(guides), kbf.exterior, n[2], True, n[1].length])

    fresh_guides_json = []
    for i, gr in enumerate(guides):
        gh = gr.copy()
        gh[2] = coords_to_flat(gh[2].coords)
        fresh_guides_json.append(gh)

    used_shapes = set([s[0] for s in guides])

    def validate_shape(shp):
        return {"id": shp[0], "length": shp[1].length, "area": Polygon(shp[1].coords).area, "type": shp[2]}

    og_poly_data = [validate_shape(pn) for pn in master_shapes if pn[0] in used_shapes]

    return guides, fresh_guides_json, og_poly_data


def merge_lines(lines: List[LineString]) -> LineString:
    last = None
    first = None
    points = []
    for line in lines:

        current = line.coords[0]

        if last is None:
            points.extend(line.coords)
        else:
            if current == last:
                points.extend(line.coords[1:])
            else:
                new_end = line.coords[-1]
                if new_end == first:
                    points[:0] = line.coords[1:]
                    # line.coords = list(line.coords)[::-1]
                else:
                    pass

        last = points[-1]
        first = points[0]

    return LineString(points)


def line_s_to_list(reso, min_length=0.0) -> List:
    result = []
    if reso.type == 'MultiLineString':
        result = [line for line in reso.geoms if line.length > min_length]
    elif reso.type == 'LineString':
        if reso.length > min_length:
            result.append(reso)
    return result


def poly_s_to_list(reso, min_length=0.0) -> List:
    result = []
    if reso.type == 'MultiPolygon':
        result = [line for line in reso.geoms if line.length > min_length]
    elif reso.type == 'Polygon':
        if reso.length > min_length:
            result.append(reso)
    return result


def get_buffered(poly, d=0.025, k=0.0075):
    ps = poly.simplify(k)
    return ps.buffer(d).simplify(k)


def get_cuts(mask, pos_poly, s_line, config, plotting=None, cut_guide=False):
    min_length = 0.0  # 0.025
    this_cut_guide = None
    if cut_guide:
        s_buf = LinearRing(get_buffered(s_line, k=config["guide_buffer"]).exterior)
        reso = s_buf.difference(pos_poly)
        buffered = line_s_to_list(reso, min_length)

        # if plotting:
        #     for line in buffered:
        #         poly_color = (0.5 * np.random.random_sample(3)).tolist()
        #         plotting.plot(*line.coords.xy, color=poly_color, linewidth=4)

        merged_guide = merge_lines(buffered)
        #print(reso.type, merged_guide.type)

        if merged_guide:
            this_cut_guide = merged_guide.intersection(mask)
            if this_cut_guide.type == 'MultiLineString':
                this_cut_guide = merge_lines(line_s_to_list(this_cut_guide))

    reso = s_line.intersection(mask)
    this_cut_line = line_s_to_list(reso, min_length)

    return this_cut_line, this_cut_guide


def simplify_poly(source_poly, config) -> List:
    simplifications = config['criteria']['simplify']
    area_limits = config['criteria']['areas']
    batch = []
    for b in range(4):
        simp_poly = source_poly.simplify(simplifications[b])
        qualify = Polygon(simp_poly)
        if qualify.area >= area_limits[b]:
            batch.append(simp_poly)
        else:
            batch.append(None)

    batch.append(source_poly)
    return batch


def parse_map_shapes(master_poly, test_bounds, config, plotting, first_pass=False):
    all_shapes = []
    guides_index = 0
    l_bounds = test_bounds.buffer(0.025, cap_style=3)
    m_bounds = test_bounds.buffer(0.05, cap_style=3)
    test_master_poly = master_poly.intersection(m_bounds)

    # if config['rect_bounded']:
    #     test_master_poly = m_bounds.difference(test_master_poly)

    # print('test_master_poly', test_master_poly.type)

    # cfe = line_s_to_list(test_master_poly)
    # les = merge_lines(cfe)
    # plotting.plot(*les.coords.xy, color="black", linewidth=2)
    # #return

    poly_set = poly_s_to_list(test_master_poly)

    minimum_outline_length_for_guide = 0.5

    for pn, poly in enumerate(poly_set):

        poly = make_valid(poly)
        show_progress('Sectors', pn, len(poly_set)+1)
        poly_color = (0.5 * np.random.random_sample(3)).tolist()
        poly_color.append(0.8)

        this_shape = {
            'id': pn,
            'bounded': True,
            'geom': [],
            'lines': [],
            'guides': []
        }

        if len(poly.interiors):
            print("\nIN interiors BOUNDS", pn)
            outline = poly.interiors[0]

            poly_outline_levels = simplify_poly(poly, config)
            outline_levels = simplify_poly(outline, config)

            this_shape['geom'] = poly_outline_levels
            this_shape['bounded'] = True
            buffered_poly = get_buffered(poly)
            these_guides = []

            for n in buffered_poly.interiors:
                spoly = Polygon(n)
                if spoly.area > 0.1:
                    this_guide = {
                        'gid': guides_index,
                        'geom': n,
                        'leng': n.length,
                        'poly': pn,
                        'closed': True
                    }
                    these_guides.append(this_guide)
                    guides_index += 1

            these_guides.sort(key=lambda x: x['leng'], reverse=True)

            this_shape['guides'] = [n for n in these_guides]

            for rn in range(5):
                if outline_levels[rn]:
                    this_line = {
                        'geom': outline_levels[rn],
                        'shape': pn,
                        'lev': rn,
                        'guides': [n['gid'] for n in these_guides],
                        'closed': True
                    }

                    this_shape['lines'].append(this_line)

        else:
            # print("\nIN BOUNDS", pn)
            outline = poly.exterior

            # print(outline)

            outline_levels = simplify_poly(outline, config)  # all LinearRings
            # print(outline_levels)
            poly_levels = simplify_poly(poly, config)
            this_shape['geom'] = poly_levels
            this_shape['bounded'] = True

            if test_bounds.contains(outline):
                guide_id = None
                if outline.length > minimum_outline_length_for_guide and 'guides' in config['includes']:
                    buffered_poly = get_buffered(poly, k=config["guide_buffer"]).exterior

                    this_guide = {
                        'gid': guides_index,
                        'geom': buffered_poly,
                        'leng': buffered_poly.length,
                        'shape': pn,
                        'closed': True
                    }
                    guide_id = this_guide['gid']
                    this_shape['guides'].append(this_guide)
                    guides_index += 1

                for rn in range(5):
                    if outline_levels[rn]:
                        this_line = {
                            'geom': outline_levels[rn],
                            'guides': [guide_id],
                            'lev': rn,
                            'shape': pn,
                            'closed': True
                        }
                        this_shape['lines'].append(this_line)

            else:
                this_shape['bounded'] = False
                for rn in range(5):
                    if outline_levels[rn]:
                        part = outline_levels[rn].intersection(l_bounds)
                        parts = line_s_to_list(part)

                        guides_or_not = rn == 4 and 'guides' in config['includes']
                        for ln, line in enumerate(parts):
                            v_lines, v_guide = get_cuts(test_bounds, poly, line, config, plotting, cut_guide=guides_or_not)

                            if v_lines and len(v_lines) and v_guide and v_guide.length > minimum_outline_length_for_guide:
                                this_guide = {
                                    'sid': ln,
                                    'gid': guides_index,
                                    'geom': v_guide,
                                    'leng': v_guide.length,
                                    'shape': pn,
                                    'closed': False
                                }
                                this_shape['guides'].append(this_guide)
                                guides_index += 1

                            this_shape['guides'].sort(key=lambda x: x['leng'], reverse=True)

                            if v_lines:
                                for v_line in v_lines:
                                    this_line = {
                                        'sid': ln,
                                        'lev': rn,
                                        'geom': v_line,
                                        'guides': [],
                                        'shape': pn,
                                        'closed': False
                                    }
                                    this_shape['lines'].append(this_line)

                for line in this_shape['lines']:
                    line['guides'] = [g['gid'] for g in this_shape['guides'] if line['sid'] == g['sid']]

        all_shapes.append(this_shape)

    return all_shapes


def OLD_get_map_shapes(master_poly, test_bounds, config, plotting, first_pass=False):
    print(master_poly.type)


    all_lines = []
    all_guides = []
    all_shapes = []

    l_bounds = test_bounds.buffer(0.025, cap_style=3)
    m_bounds = test_bounds.buffer(0.1, cap_style=3)
    test_master_poly = master_poly.intersection(m_bounds)

    if config['rect_bounded']:
        test_master_poly = m_bounds.difference(test_master_poly)

    print(test_master_poly.type)

    # #// in the case of the corsica contour there are some randy extra polys.
    # #// presence of interior indicated enclosing parameter.

    buf_depth = 0.0125
    min_guide_area = 0.01
    minimum_outline_length_for_guide = 0.1

    if test_master_poly.type == 'MultiPolygon':
        # #//each positive polygon here
        for pn, poly in enumerate(test_master_poly.geoms):
            show_progress('Sectors', pn, len(test_master_poly.geoms))
            poly_color = (0.5 * np.random.random_sample(3) ).tolist()
            poly_color.append(0.8)

            this_shape = {
                'id': pn,
                'bounded': True,
                'geom': poly,
                'lines': [],
                'guides': [],
                'levels': []
            }

            """
            SHAPE {LINES, GUIDES}
            
            when dealing with an enclosed outer shape, like the whole med, 
            the line has to be collected from inside the thing.
            if that main contour...
            do we one pass at max resolution for the guides, 
            they are indexed per shape/path
            
            make all then associated subsequent iterations.
            line index stays the same.
            
            make first pass, make lower-res passes
            """

            if len(poly.interiors):
                # this_line = {
                #     'id': pn,
                #     'geom': poly.interiors[0],
                #     'guides': []
                # }
                outline = poly.interiors[0]
                poly_outline_levels = simplify_poly(poly, config)  # all LinearRings
                outline_levels = simplify_poly(outline, config)  # all LinearRings

                this_shape['geom'] = poly_outline_levels
                this_shape['bounded'] = True

                buffered_poly = get_buffered(poly)
                these_guides = []
                start_index = len(all_guides)
                for n in buffered_poly.interiors:
                    spoly = Polygon(n)
                    if spoly.area > 0.1:
                        this_guide = {
                            'id': len(these_guides)+start_index,
                            'geom': n,
                            'leng': n.length,
                            'poly': pn
                        }
                        these_guides.append(this_guide)

                #         plotting.plot(*n.coords.xy, color=(1.0, 0.0, 1.0, 0.9), linewidth=1)
                #
                # for ote in outline_levels:
                #     if ote:
                #         plotting.plot(*ote.coords.xy, color=poly_color, linewidth=1)

                # plotting.plot(*this_line['geom'].coords.xy, color=(0.0, 0.0, 0.0, 0.9), linewidth=1)

                these_guides.sort(key=lambda x: x['leng'], reverse=True)

                this_shape['guides'] = [n for n in these_guides]  #.append(this_guide)

                this_line = {
                    'id': len(this_shape['lines']),
                    'geom': outline_levels,
                    'shape': pn
                }

                this_shape['lines'].append(this_line)
                #
                #
                # guide_indices = [n['id'] for n in these_guides]
                # for g in these_guides:
                #     all_guides.append(g)
                #
                # this_line['guides'] = guide_indices
                # all_lines.append(this_line)

                # print(all_lines)
                # print(all_guides)

            else:

                outline = poly.exterior
                outline_levels = simplify_poly(outline, config)  # all LinearRings
                this_shape['geom'] = outline_levels
                this_shape['bounded'] = True

                if test_bounds.contains(outline):
                    # #// in bounds so shape (geom)(5) and guide(1)
                    # for ote in outline_levels:
                    #     if ote:
                    #         plotting.plot(*ote.coords.xy, color=poly_color, linewidth=1)

                    guide_id = None
                    if outline.length > minimum_outline_length_for_guide:
                        buffered_poly = get_buffered(poly).exterior
                        this_guide = {
                            'id': len(this_shape['guides']),
                            'geom': buffered_poly,
                            'leng': buffered_poly.length,
                            'shape': pn
                        }
                        guide_id = this_guide['id']
                        this_shape['guides'].append(this_guide)

                        # plotting.plot(*this_guide['geom'].coords.xy, color=poly_color, linewidth=1)

                    this_line = {
                        'id': len(this_shape['lines']),
                        'geom': outline_levels,
                        'guide': [guide_id],
                        'shape': pn
                    }

                    this_shape['lines'].append(this_line)

                else:
                    # #// out of bounds so shape (geom)(5) and guide(x)
                    this_shape['bounded'] = False
                    for rn in range(5):
                        if outline_levels[rn]:
                            part = outline_levels[rn].intersection(l_bounds)
                            parts = line_s_to_list(part)

                            for ln, line in enumerate(parts):
                                # #//get cuts are only useful for guides creation.
                                v_lines, v_guide = get_cuts(test_bounds, poly, line, cut_guide=rn == 4)

                                if v_lines and len(v_lines) and v_guide:
                                    this_guide = {
                                        'id': len(this_shape['guides']),
                                        'sid': ln,
                                        'geom': v_guide,
                                        'leng': v_guide.length,
                                        'shape': pn
                                    }
                                    # has_guide = this_guide['id']
                                    this_shape['guides'].append(this_guide)
                                    #plotting.plot(*v_guide.coords.xy, color=poly_color, linewidth=4)

                                if v_lines:
                                    for v_line in v_lines:
                                        this_line = {
                                            'id': len(this_shape['lines']),
                                            'sid': ln,
                                            'lev': rn,
                                            'geom': v_line,
                                            'poly': pn
                                        }
                                        this_shape['lines'].append(this_line)
                                        #plotting.plot(*v_line.coords.xy, color='black', linewidth=2)



            #plot_polygon(plotting, poly, color=poly_color)
            all_shapes.append(this_shape)

            if pn == 0:
                # if not this_shape['bounded'] and len(this_shape['guides']) > 1:

                out = json.dumps(this_shape, indent=2, cls=JsonSafeEncoder)
                print(out)

            # print(this_shape['id'])
            # print(this_shape['geom'])
            # print(this_shape['lines'])

            # for k, m_shape in this_shape.items():
            #
            #
            #
            #
            #     if type(v) == list:
            #         p = len(v)
            #     elif type(v) == int:
            #         p = v
            #     else:
            #         p = type(v)
            #     print(k, p)

    # print(all_lines)
    # print(all_guides)


def new_data_from_sector(src_part, src, sector_tuple, sector_width, offset=(0, 0), mod=1) -> np.array:
    #// evaluate orientation of lats first: consistent issue depending on dataset..assume top down
    # lat, lon
    # time, lat, lon
    # time, depth, lat, lon
    dimensions = len(src_part.shape)
    ye = np.where((src['lats'] <= sector_tuple[1]) & (src['lats'] > (sector_tuple[1] - sector_width)))[0]  # indices
    xe = np.where((src['lons'] >= sector_tuple[0]) & (src['lons'] < (sector_tuple[0] + sector_width)))[0]  # indices

    if dimensions == 2:
        result = src['data'][ye[offset[1]:, None], xe+offset[0]]
    elif dimensions == 3:
        result = src_part[:, ye[offset[1]:, None], xe+offset[0]]
    elif dimensions == 4:
        result = src['data'][:, :, ye[offset[1]:, None], xe+offset[0]]
    else:
        result = None

    return result
