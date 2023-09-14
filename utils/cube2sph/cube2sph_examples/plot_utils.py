import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
from xyz2latlond import xyz2latlond
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
def get_taper(t, t_start, t_end, order=1):
    tw = np.ones_like(t, np.float)
    ind = (t >= t_start)
    tw[ind] = np.power(np.cos((t[ind] - t_start) / (t_end - t_start) * np.pi / 2.0), order)
    ind = (t >= t_end)
    tw[ind] = 0.0
    return tw

def get_taper_2side(t, t1, t2, t3, t4, order=1):
    tw = np.ones_like(t, np.float)
    tw[t>t4] = 0.0
    tw[t<t1] = 0.0
    ind = np.logical_and((t > t2), (t < t3))
    tw[ind] = 1.0
    ind = np.logical_and((t >= t3), (t <= t4))
    tw[ind] = np.power(np.cos((t[ind] - t3) / (t4 - t3) * np.pi / 2.0), order)
    ind = np.logical_and((t >= t1), (t <= t2))
    tw[ind] = np.power(np.cos((t[ind] - t2) / (t2 - t1) * np.pi / 2.0), order)
    return tw

def split_line(line):
  # split line with continuous multiple spaces
  # return a list of strings
  return [_ for _ in line.strip().split(' ') if _ != '']

def get_dist_deg_az(coord1, coord2, is_geo_coord=False, ASSUME_PERFECT_SPHERE=False):
  # if is_geo_coord is set to True, coord=(lat, lon),
  # otherwise, coord = (x, y, z)
  if is_geo_coord:
    lat1 = coord1[0]
    lon1 = coord1[1]
    lat2 = coord2[0]
    lon2 = coord2[1]
  else:
    x1 = coord1[0]
    y1 = coord1[1]
    z1 = coord1[2]
    x2 = coord2[0]
    y2 = coord2[1]
    z2 = coord2[2]
    lat1, lon1, d1 = xyz2latlond(x1, y1, z1, ASSUME_PERFECT_SPHERE)
    lat2, lon2, d2 = xyz2latlond(x2, y2, z2, ASSUME_PERFECT_SPHERE)
  if ASSUME_PERFECT_SPHERE:
    dist, az, baz = gps2dist_azimuth(lat1, lon1, lat2, lon2, 6371000.0, 0.0)
  else:
    dist, az, baz = gps2dist_azimuth(lat1, lon1, lat2, lon2)
  dist = kilometers2degrees(dist / 1000.0)
  return dist, az, baz



class WaveformSection:
  def __init__(self):
    self.waveforms = None
    self.t = None
    self.dist_list = None
    self.az_list = None
    self.baz_list = None
    self.windows = None

  def __sub__(self, other):
    sub = other.interpolate(self.t)
    for i_wave in range(len(sub.waveforms)):
      sub.waveforms[i_wave] = self.waveforms[i_wave] - sub.waveforms[i_wave]
    return sub

  def difference_ratio(self, other):
    ratio_list = []
    sub = other.interpolate(self.t)
    for i_wave in range(len(sub.waveforms)):
      ind = np.logical_and((self.t >= self.windows[i_wave][0]), (self.t <= self.windows[i_wave][1]))
      diff = self.waveforms[i_wave][ind] - sub.waveforms[i_wave][ind]
      ratio = np.sum(diff * diff) / (
              np.sqrt(np.sum(self.waveforms[i_wave][ind] * self.waveforms[i_wave][ind]))*
              np.sqrt(np.sum(sub.waveforms[i_wave][ind] * sub.waveforms[i_wave][ind])))
      ratio_list.append(ratio)
    return ratio_list

  def time_shift(self, other):
    shift_list = []
    sub = other.interpolate(self.t)
    for i_wave in range(len(sub.waveforms)):
      ind = np.logical_and((self.t >= self.windows[i_wave][0]), (self.t <= self.windows[i_wave][1]))
      cc = signal.correlate(self.waveforms[i_wave][ind], sub.waveforms[i_wave][ind])
      dt = (np.argmax(cc) - (len(self.waveforms[i_wave][ind]) - 1)) * (self.t[1] - self.t[0])
      shift_list.append(dt)
    return shift_list

  def get_diff_gram(self, other, tmin, dt, twin, amp_threshold=0.01):
    tshift_gram_list = []
    cc_gram_list = []
    tc_gram_list = []
    sub = other.interpolate(self.t)
    for i_wave in range(len(sub.waveforms)):
      tc = tmin
      tshift_gram = []
      cc_gram = []
      tc_gram = []
      ind_wave = np.logical_and((self.t >= self.windows[i_wave][0]), (self.t <= self.windows[i_wave][1]))
      max_amp = np.amax(np.abs(self.waveforms[i_wave][ind_wave]))
      while (tc+twin/2.0 < self.windows[i_wave][1]):
        tc_gram.append(tc)
        ind_win = np.logical_and((self.t >= tc-twin/2.0), (self.t <= tc+twin/2.0))
        wave1 = self.waveforms[i_wave][ind_win]
        wave2 = sub.waveforms[i_wave][ind_win]
        if (np.amax(np.abs(wave1)) < max_amp * amp_threshold or np.amax(np.abs(wave2)) < max_amp * amp_threshold):
          tshift_gram.append(None)
          cc_gram.append(None)
          tc = tc + dt
          continue
        energy1 = np.sqrt(np.sum(wave1*wave1))
        energy2 = np.sqrt(np.sum(wave2*wave2))
        cc = signal.correlate(wave1, wave2) / energy1 / energy2
        cc_gram.append(np.amax(cc))
        tshift_gram.append((np.argmax(cc) - (len(wave1) - 1)) * (self.t[1] - self.t[0]))
        tc = tc + dt
      tshift_gram_list.append(tshift_gram)
      cc_gram_list.append(cc_gram)
      tc_gram_list.append(tc_gram)
    return tc_gram_list, tshift_gram_list, cc_gram_list



  def scale(self, fac):
    for i_wave in range(len(self.waveforms)):
      self.waveforms[i_wave] = self.waveforms[i_wave] * fac

  def set_windows_default(self):
    self.windows = []
    for i_wave in range(len(self.waveforms)):
      self.windows.append((self.t.min(), self.t.max()))

  def set_windows(self, vmin=None, vmax=None, tmin=None, tmax=None, in_degree=True):
    if tmin is None: tmin = self.t.min() 
    if tmax is None: tmax = self.t.max()
    if not in_degree:
      # velocity in km/s, not degree/s
      if vmin is not None: vmin = kilometers2degrees(vmin)
      if vmax is not None: vmax = kilometers2degrees(vmax)
    self.windows = []
    for i_wave in range(len(self.waveforms)):
      if (vmax is None):
        t_min = max(tmin, self.t.min())
      else:
        t_min = max(self.dist_list[i_wave] / vmax, tmin, self.t.min())
      if (vmin is None):
        t_max = min(tmax, self.t.max())
      else:
        t_max = min(self.dist_list[i_wave] / vmin, tmax, self.t.max())
      if (t_max <= (t_min + self.t[1] - self.t[0])):
        sys.exit(f"({t_min}, {t_max}) is not a valid window\n")
      self.windows.append((t_min, t_max))
      

  def interpolate(self, t_new):
    """
    not in-place
    """
    sec_new = WaveformSection()
    sec_new.dist_list = self.dist_list
    sec_new.az_list = self.az_list
    sec_new.baz_list = self.baz_list
    sec_new.t = np.copy(t_new)
    sec_new.waveforms = []
    for wave in self.waveforms:
      f = interp1d(self.t, wave, fill_value="extrapolate")
      sec_new.waveforms.append(f(t_new))
    sec_new.windows = []
    for win in self.windows:
      sec_new.windows.append(win)
    return sec_new

  def taper(self, taper_range, verbose=1):
    if (not isinstance(taper_range, tuple)):
      sys.exit(f"{taper_range} is not a tuple\n")
    if (len(taper_range)==2):
      if (verbose==1): print(f'getting 2-point taper {taper_range}\n')
      tw = get_taper(self.t, taper_range[0], taper_range[1])
    elif (len(taper_range)==4):
      if (verbose==1): print(f'getting 4-point taper {taper_range}\n')
      tw = get_taper_2side(self.t, taper_range[0], taper_range[1], 
                                   taper_range[2], taper_range[3])
    else:
      sys.exit(f"incorrect taper {taper_range} \n")
    for i_wave in range(len(self.waveforms)):
      wave = self.waveforms[i_wave]
      self.waveforms[i_wave] = wave * tw
  
  def filter(self, freq_range, verbose=1):
    if (not isinstance(freq_range, tuple)):
      sys.exit(f"{freq_range} is not a tuple\n")
    if (len(freq_range)!=2):
      sys.exit(f"incorrect frequency range {freq_range}\n")
    if (verbose==1): print(f'apply filter {freq_range}\n')
    s = self.t[1] - self.t[0]
    sos = signal.butter(4, [freq_range[0] * s * 2, freq_range[1] * s * 2], 'bandpass', output='sos')
    for i_wave in range(len(self.waveforms)):
      wave = self.waveforms[i_wave]
      self.waveforms[i_wave] = signal.sosfiltfilt(sos, wave, padtype=None)

  def plot_waveforms(self, *args, fig_num=0, time_range=None, normalize=1.0, 
                     norm_with=None, label=None, **kwargs):
    fig = plt.figure(num=fig_num) # pull out the existed figure
    #if time_range is None: time_range = (self.t[0], self.t[-1])
    #ind = np.logical_and((self.t >= time_range[0]), (self.t <= time_range[1]))
    for i_wave in range(len(self.waveforms)):
      if (time_range is None):
        ind = np.logical_and((self.t >= self.windows[i_wave][0]), (self.t <= self.windows[i_wave][1]))
      else:
        ind = np.logical_and((self.t >= time_range[0]), (self.t <= time_range[1]))
      wave = self.waveforms[i_wave][ind]
      t_ind = self.t[ind]
      if (norm_with is None):
        norm_val = np.amax(wave)
      else:
        if (time_range is None):
          ind_norm = np.logical_and((norm_with.t >= norm_with.windows[i_wave][0]), (norm_with.t <= norm_with.windows[i_wave][1]))
        else:
          ind_norm = np.logical_and((norm_with.t >= time_range[0]), (norm_with.t <= time_range[1]))
        norm_val = np.amax(norm_with.waveforms[i_wave][ind_norm])
      wave_norm = wave / norm_val * normalize + self.dist_list[i_wave]
      if ((i_wave == 0) and (label is not None)):
        plt.plot(t_ind, wave_norm, *args, label=label, **kwargs)
      else:
        plt.plot(t_ind, wave_norm, *args, **kwargs)

  def plot_diff_gram(self, tc_gram_list, diff_gram_list, *args, fig_num=0, 
                     normalize=1.0, label=None, **kwargs):
    fig = plt.figure(num=fig_num) # pull out the existed figure
    add_label = True
    for i_wave in range(len(self.waveforms)):
      for it in range(len(tc_gram_list[i_wave])):
        if (diff_gram_list[i_wave][it] is not None):
          if add_label and (label is not None):
            plt.plot(tc_gram_list[i_wave][it], diff_gram_list[i_wave][it] * normalize + self.dist_list[i_wave], *args, label=label, **kwargs)
            add_label = False
          else:
            plt.plot(tc_gram_list[i_wave][it], diff_gram_list[i_wave][it] * normalize + self.dist_list[i_wave], *args, **kwargs)



  def get_waveforms_from_axisem(self,
                                fn_list, 
                                t=None,
                                tstart=None,
                                dt=None,
                                nt=None,
                                fn_source=None, 
                                fn_stations=None,
                                verbose=1):
  
    if (verbose==1): print('plotting waveforms from AxiSEM results\n')
    n_stations = len(fn_list)
    if (verbose==1): print(f'there are {n_stations} waveforms\n')
    if t is None:
      if ((tstart is None) or (dt is None) or (nt is None)):
        wave = np.loadtxt(fn_list[0])
        if (tstart is None):
          if (verbose==1): print(f'using tstart from {fn_list[0]}')
          tstart = wave[0,0]
        if (dt is None): 
          if (verbose==1): print(f'using dt from {fn_list[0]}')
          dt = wave[1,0] - wave[0,0]
        if (nt is None): 
          if (verbose==1): print(f'using nt from {fn_list[0]}')
          nt = wave.shape[0]
      if (verbose==1): print(f'simulation time: tstart={tstart}, dt={dt}, nt={nt}\n')
      try:
        t = np.arange(0, nt, dtype=float) *dt + tstart
      except Exception as e:
        sys.exit(f"{e}\n")
    self.t = t

    if (verbose==1): print(f'reading source file {fn_source}\n')
    try:
      with open(fn_source, 'r') as f_source:
        lines = f_source.readlines()
      for line in lines:
        if (line.lstrip().startswith('SOURCE_LAT') or line.lstrip().startswith('lat')):
          source_lat = float(split_line(line)[1])
        if (line.lstrip().startswith('SOURCE_LON') or line.lstrip().startswith('lon')):
          source_lon = float(split_line(line)[1])
      dist_list = [None] * n_stations
      az_list = [None] * n_stations
      baz_list = [None] * n_stations
      with open(fn_stations, 'r') as f_stations:
        lines = f_stations.readlines()
      for line in lines:
        line_segs = split_line(line)
        nt = line_segs[0]
        sta = line_segs[1]
        sta_lat = float(line_segs[2])
        sta_lon = float(line_segs[3])
        for i_fn in range(len(fn_list)):
          fn = os.path.split(fn_list[i_fn])[-1]
          if (fn.startswith(f'{nt}_{sta}') or fn.startswith(f'{sta}_{nt}')):
            if dist_list[i_fn] is not None:
              raise(ValueError('file list conflict'))
            dist_list[i_fn], az_list[i_fn], baz_list[i_fn] = get_dist_deg_az(
                                    (source_lat, source_lon), 
                                    (sta_lat, sta_lon), 
                                    is_geo_coord = True,
                                    ASSUME_PERFECT_SPHERE=True)
            continue
      for i_fn in range(len(fn_list)):
        if dist_list[i_fn] is None:
          raise(ValueError(f'station does not exist for {fn_list[i_fn]}'))
    except Exception as e:
      sys.exit(f"{e}\n")
    self.dist_list = dist_list
    self.az_list = az_list
    self.baz_list = baz_list

    waveforms = []
    for i_fn in range(len(fn_list)):
      try:
        if (verbose==1): print(f'reading waveform file {fn_list[i_fn]}\n')
        wave = np.loadtxt(fn_list[i_fn])
        if (wave.shape[0] != len(t)):
          raise (ValueError(f'file length inconsistant: {fn_list[i_fn]}'))
        waveforms.append(-wave)
      except Exception as e:
        sys.exit(f"{e}\n")
    self.waveforms = waveforms
    self.set_windows_default()

  def get_waveforms_from_specfem(self,
                                fn_list,
                                t=None,
                                tstart=None,
                                dt=None,
                                nt=None,
                                fn_source=None,
                                fn_stations=None,
                                is_cartesian=True,
                                ellipticity=True,
                                delimiter='.',
                                verbose=1):

    if (verbose==1): print('plotting waveforms from SPECFEM results\n')
    n_stations = len(fn_list)
    if (verbose==1): print(f'there are {n_stations} waveforms\n')
    if t is None:
      if ((tstart is None) or (dt is None) or (nt is None)):
        wave = np.loadtxt(fn_list[0])
        if (tstart is None):
          if (verbose==1): print(f'using tstart from {fn_list[0]}')
          tstart = wave[0,0]
        if (dt is None): 
          if (verbose==1): print(f'using dt from {fn_list[0]}')
          dt = wave[1,0] - wave[0,0]
        if (nt is None): 
          if (verbose==1): print(f'using nt from {fn_list[0]}')
          nt = wave.shape[0]
      if (verbose==1): print(f'simulation time: tstart={tstart}, dt={dt}, nt={nt}\n')
      try:
        t = np.arange(0, nt, dtype=float) *dt + tstart
      except Exception as e:
        sys.exit(f"{e}\n")
    self.t = t

    if (verbose==1): print(f'reading source file {fn_source}\n')
    try:
      with open(fn_source, 'r') as f_source:
        lines = f_source.readlines()
      for line in lines:
        if line.lstrip().startswith('lat'):
          source_lat = float(split_line(line)[1])
        if line.lstrip().startswith('lon'):
          source_lon = float(split_line(line)[1])
        if line.lstrip().startswith('depth'):
          source_depth = float(split_line(line)[1])
      dist_list = [None] * n_stations
      az_list = [None] * n_stations
      baz_list = [None] * n_stations
      with open(fn_stations, 'r') as f_stations:
        lines = f_stations.readlines()
      for line in lines:
        line_segs = split_line(line)
        nt = line_segs[0]
        sta = line_segs[1]
        sta_lat = float(line_segs[2])
        sta_lon = float(line_segs[3])
        sta_depth = float(line_segs[5])
        for i_fn in range(len(fn_list)):
          fn = os.path.split(fn_list[i_fn])[-1]
          if (fn.startswith(f'{nt}{delimiter}{sta}') or fn.startswith(f'{sta}{delimiter}{nt}')):
            if dist_list[i_fn] is not None:
              raise(ValueError('file list conflict'))
            if is_cartesian:
              dist_list[i_fn], az_list[i_fn], baz_list[i_fn] = get_dist_deg_az(
                                    (source_lon, source_lat, source_depth),
                                    (sta_lon, sta_lat, sta_depth),
                                    ASSUME_PERFECT_SPHERE = (not ellipticity), 
                                    is_geo_coord = False)
            else:
              dist_list[i_fn], az_list[i_fn], baz_list[i_fn] = get_dist_deg_az(
                                    (source_lat, source_lon),
                                    (sta_lat, sta_lon),
                                    ASSUME_PERFECT_SPHERE = (not ellipticity),
                                    is_geo_coord = True)
            continue
      for i_fn in range(len(fn_list)):
        if dist_list[i_fn] is None:
          raise(ValueError(f'station does not exist for {fn_list[i_fn]}'))
    except Exception as e:
      sys.exit(f"{e}\n")
    self.dist_list = dist_list
    self.az_list = az_list
    self.baz_list = baz_list

    waveforms = []
    for i_fn in range(len(fn_list)):
      try:
        if (verbose==1): print(f'reading waveform file {fn_list[i_fn]}\n')
        wave = np.loadtxt(fn_list[i_fn])
        if (wave.shape[0] != len(t)):
          raise (ValueError(f'file length inconsistant: {fn_list[i_fn]}'))
        waveforms.append(wave[:,1])
      except Exception as e:
        sys.exit(f"{e}\n")
    self.waveforms = waveforms
    self.set_windows_default()


def rotate_ne_rt(n, e, a):
  a = a / 180.0 * np.pi
  r = - e * np.sin(a) - n * np.cos(a)
  t = - e * np.cos(a) + n * np.sin(a)
  return r, t

def rotate_section_ne2rt_receiver(sec_n, sec_e):
  assert len(sec_n.waveforms) == len(sec_e.waveforms), "rotate error: inconsistent number of waveforms"
  sec_r = WaveformSection()
  sec_t = WaveformSection()
  sec_r.t = sec_n.t
  sec_r.dist_list = sec_n.dist_list
  sec_r.az_list = sec_n.az_list
  sec_r.baz_list = sec_n.baz_list
  sec_r.waveforms = []
  sec_t.t = sec_n.t
  sec_t.dist_list = sec_n.dist_list
  sec_t.az_list = sec_n.az_list
  sec_t.baz_list = sec_n.baz_list
  sec_t.waveforms = []
  for i_wave in range(len(sec_n.waveforms)):
    assert abs(sec_n.baz_list[i_wave] - sec_e.baz_list[i_wave]) < 1.0e-3, "rotate error: inconsistent backazimuth"
    n = sec_n.waveforms[i_wave]
    e = sec_e.waveforms[i_wave]
    ba = sec_n.baz_list[i_wave]
    r, t = rotate_ne_rt(n, e, ba)
    sec_r.waveforms.append(r)
    sec_t.waveforms.append(t)
  sec_r.set_windows_default()
  sec_t.set_windows_default()
  return sec_r, sec_t

def rotate_section_ne2rt_source(sec_n, sec_e):
  assert len(sec_n.waveforms) == len(sec_e.waveforms), "rotate error: inconsistent number of waveforms"
  sec_r = WaveformSection()
  sec_t = WaveformSection()
  sec_r.t = sec_n.t
  sec_r.dist_list = sec_n.dist_list
  sec_r.az_list = sec_n.az_list
  sec_r.baz_list = sec_n.baz_list
  sec_r.waveforms = []
  sec_t.t = sec_n.t
  sec_t.dist_list = sec_n.dist_list
  sec_t.az_list = sec_n.az_list
  sec_t.baz_list = sec_n.baz_list
  sec_t.waveforms = []
  for i_wave in range(len(sec_n.waveforms)):
    assert abs(sec_n.az_list[i_wave] - sec_e.az_list[i_wave]) < 1.0e-3, "rotate error: inconsistent azimuth"
    n = sec_n.waveforms[i_wave]
    e = sec_e.waveforms[i_wave]
    az = sec_n.az_list[i_wave] 
    r, t = rotate_ne_rt(n, e, az+180.0)
    sec_r.waveforms.append(r)
    sec_t.waveforms.append(t)
  sec_r.set_windows_default()
  sec_t.set_windows_default()
  return sec_r, sec_t
