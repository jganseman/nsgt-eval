#!/usr/bin/python

import sys,time

def datesuffix():
    # Create a time suffix.
    t=time.localtime()
    
    date=`t[2]`
    if len(date)==1:
        date='0'+date
        
    month=`t[1]`
    if len(month)==1:
        month='0'+month

    year=`t[0]`[2:]

    s = date+month+year

    return s

if 'datesuffix' in sys.argv[0]:
    print datesuffix()
