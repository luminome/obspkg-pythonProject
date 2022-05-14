import json
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon, Point, box, LineString, LinearRing
from shapely.ops import nearest_points


dpaths = [
    ['/Users/sac/Downloads/wudi_data.json', 'wudi'],
    ['/Users/sac/Downloads/guideo_data-2.json', 'guides']
]

group = {
    'guides': [],
    'wudi': [],
}

guides = []

fig, ax = plt.subplots()

for e in dpaths:
    with open(e[0]) as f:
        group[e[1]] = json.load(f)

for guide in group['guides']:
    data_part = np.array(guide['data'])[:, 1:3]


    print(data_part.shape)
    M = LineString(data_part)
    guides.append(M)
    plt.plot(*M.coords.xy)


wudi_part = np.array(group['wudi'])[:, 1:3]
M = LineString(wudi_part)
plt.plot(*M.coords.xy)
plt.scatter(wudi_part[:, 0], wudi_part[:, 1], s=3)

sta = []
guide = guides[0]

for p in M.coords:
    # pt = Point(p)
    # p2 = nearest_points(guide, pt)[0]
    # print(p2)  # POINT (5 7)
    # sta.append((p2.x, p2.y))

    pt = Point(p)
    p2 = guide.interpolate(guide.project(pt))
    sta.append((p2.x, p2.y))
#
#
    d = pt.distance(p2)
    if d < 0.1:
        N = np.array([[pt.x, pt.y], [p2.x, p2.y]])
        plt.plot(N[:, 0], N[:, 1], marker='o')
#
fta = np.array(sta)
print(fta)

plt.scatter(fta[:, 0], fta[:, 1], s=4)



plt.grid()
plt.tight_layout()
plt.axis('equal')
plt.show()





# out = json.dumps(group, indent=2)
# print(out)
exit()
#
#
#     group[e[1]] = np.array(group[e[1]])
#
#     print(group[e[1]].shape)
#
#     ft = group[e[1]][:, 1:3]
#     print(ft.shape)
#
#     M = LineString(ft)
#
#
#
#     plt.plot(*M.coords.xy)
#
#     # # plt.plot(M[:, 0], M[:, 1])  # , s=4.0, marker='o')
#     #
#     # plt.scatter(my_data['lon'], my_data['lat'], s=1)
#     #
#     # plt.plot(my_data['lonmid'], my_data['latmid'], marker='+')
#     #
#     # plt.plot(my_data['lon'], my_data['lat'], marker='o')
#     #
#     # out = tool.do_refine(my_data['lon'], my_data['lat'], 1, 3)
#     # plt.plot(out[0], out[1], linewidth=1, color=(0.0, 0.0, 0.0))
#
#
#
# plt.grid()
# plt.tight_layout()
# plt.axis('equal')
# plt.show()


# out = json.dumps(group, indent=2)
# print(out)