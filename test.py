#!/usr/bin/python
import psycopg2
import json


"""
#do not:
#DROP TABLE IF EXISTS sample_geolocs;
#instead:
#DELETE FROM sample_geolocs;



"""
#//TODO:SELECT id,shape,hole FROM sample_geolocs WHERE id % 5 = 0 LIMIT 100;

# conn = "ok"
conn = psycopg2.connect(
    host="localhost",
    database="sac")


def do_load():
    # Opening JSON file
    f = open('/Users/sac/downloads/data-24.json')  #;//./data/data-3 2.json')

    # returns JSON object as
    # a dictionary
    data = json.load(f)
    return data


def insert_point(cursor, connection, point_json):
    sql = """INSERT INTO sample_geolocs(lon,lat,shape) VALUES(%s,%s,%s) RETURNING id;"""
    try:
        variable = 'NULL'
        #hole = None if point_json['hole'] == 'null' else point_json['hole']
        prep = (point_json['v']['lon'], point_json['v']['lat'], point_json['shape'])  #, hole,)
        #text = text.replace("nan", "null")

        cursor = connection.cursor()
        print(prep)
        # execute the INSERT statement
        cursor.execute(sql, prep)
        # get the generated id back
        the_id = cursor.fetchone()[0]
        print('the_id', the_id)
        # commit the changes to the database
        connection.commit()
        # close communication with the database
        #cursor.close()
    except (Exception, psycopg2.DatabaseError) as error:
        print('error', error)
    # finally:
    #     if connection is not None:
    #         connection.close()


def get_one(connection, res):
    sql = """SELECT lon, lat, shape FROM sample_geolocs WHERE mod(id,%s) = 0 ORDER BY shape ASC;"""
    try:
        cursor = connection.cursor()
        cursor.execute(sql, res)
        data = cursor.fetchall()
        return data
    except (Exception, psycopg2.DatabaseError) as error:
        print('error', error)


def deliver(data, name):
    filename = f'data/data-{name}-test.json'
    with open(filename, 'w') as fp:
        json.dump(data, fp, indent=2)


if __name__ == '__main__':
    print(conn)
    # create a cursor
    cur = conn.cursor()
    # execute a statement
    print('PostgreSQL database version:')
    cur.execute('SELECT version()')
    # display the PostgreSQL database server version
    db_version = cur.fetchone()
    print(db_version)

    dat = get_one(conn, [10])
    #print(dat)
    whole_lex = []
    for v, i in enumerate(dat):
        print(v, i)
        whole_lex.append(i)

    deliver(whole_lex, 'sample')

    # close the communication with the PostgreSQL
    #//THIS IS THE INSERTION SCRIPT:
    #//cur.close()

    # dat = do_load()
    # for v, i in enumerate(dat):
    #     print(v, i)
    #     insert_point(cur, conn, i)

    #//TODO: 48996 points total


#//TODO: how to sectorize the delivery?
#//shapes?
