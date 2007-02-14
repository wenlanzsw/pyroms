# encoding: utf-8

# cdftime
# Copyright (c) 2006 Robert Hetland


from numpy import *
from netcdftime import utime
import netCDF4, datetime
from pylab import date2num
from pyroms import Dataset

class roms_time (ndarray):
    
    def __new__(self, ncfile, name='ocean_time', units=None, **kwargs):
        self._nc = Dataset(ncfile)
        romstime = self._nc.variables[name][:]
        if units == None:
            self._units = self._nc.variables[name].units
        else:
            self._units = units
        self._utime = utime(self._units, **kwargs)
        self.origin = self._utime.origin
        self.dates = self._utime.num2date(romstime)
        
        return romstime.view(roms_time)
    
    def nearest_index(self, dateo):
        to = self._utime.date2num(dateo)
        print 'to = ', to
        t = self._utime.date2num(self.dates)
        print 't = ', t
        return where(abs(t-to) == min(abs(t-to)))
    
    def nearest(self, dateo):
        return self.dates[self.nearest_index(dateo)]
    
    def get_seconds(self):
        self._utime.units = 'seconds'
        return self._utime.date2num(self.dates)
    
    def get_minutes(self):
        self._utime.units = 'minutes'
        return self._utime.date2num(self.dates)
    
    def get_hours(self):
        self._utime.units = 'hours'
        return self._utime.date2num(self.dates)
    
    def get_days(self):
        self._utime.units = 'seconds'
        return self._utime.date2num(self.dates)
    
    def get_jd(self):
        return date2num(self.dates)
        
    jd = property(get_jd, None, doc="Julian day, for plotting in pylab")
    seconds = property(get_seconds, None, doc="seconds")
    minutes = property(get_minutes, None, doc="minutes")
    hours = property(get_hours, None, doc="hours")
    days = property(get_days, None, doc="days")

if __name__ == '__main__':
    nc = netCDF4.Dataset('/Volumes/Scratch/near_5_250.0_his.nc')
    time = time(nc,time='ocean_time')
    print time.dates
    print time.dates[4:8]
    print time._utime.units
    print time
    print time[2:8]
    print time.seconds
    print time.days
    print time.hours
    do = datetime.datetime(1,01,01)
    print time.nearest_index(do)
    print time.nearest(do)
    
    print time[:]
    print len(time)
    
