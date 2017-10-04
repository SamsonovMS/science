"""
Write constants to txt.
"""
import c

def write_constants(**kwargs):
    """
    Create txt with constants from animation.py.
    """
    counter = 0
    log = open('constants.txt', 'a')
    for name in c.constants_names:
        if name in kwargs:
            counter += 1
            log.write(name + ' ' + str(kwargs[name]) + '\n')
        else:
            print("write_constants: no %s in kwargs" % name)
            print(kwargs)
    log.close()
    if counter != len(kwargs):
        print("%s != len(kwargs)" % counter)
        print(kwargs)