# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 13:54:23 2015

@author: Weberjm
"""
try:
    # Python 3
    from functools import reduce
except:
    pass

def numToTime(time, unit='s', outunit='time'):
    '''
    Function to convert time with a unit to another unit.
    '''
    unit = unit.lower()

    time = float(time)

    # convert time to seconds
    if unit in ['d','days']:
        time *= 60 * 60 * 24
    elif unit in ['h','hr']:
        time *= 60 * 60
    elif unit in ['m', 'min']:
        time *= 60

    if outunit=='time':
        time = reduce(lambda ll,b : divmod(ll[0],b) + ll[1:], [(time,),60,60,24])

        tlist = []
        for i, num in enumerate(time):
            if i == 3:
                tlist.append('{:.2f}'.format(num))
            elif num>0 or len(tlist)>0:
                tlist.append(int(num))

        return ':'.join([str(t) for t in tlist])

    elif outunit in ['d', 'days']:
        return time/(60.0*60.0*24.0)

    elif outunit in ['h', 'hr', 'hrs']:
        return time/(60.0*60.0)

    elif outunit in ['m', 'min', 'mins']:
        return time/(60.0)

    else:
        return time