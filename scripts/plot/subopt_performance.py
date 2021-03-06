# Copyright 2016, Eliot Courtney.
#
# This file is part of kekrna.
#
# kekrna is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# kekrna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with kekrna.
# If not, see <http://www.gnu.org/licenses/>.
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from plot.fold_accuracy import TEXT_LOC
from plot.load_data import colmap
from plot.plot_common import savefig_local, do_quantity_plot, do_quantity_log_plot, get_subplot_grid, \
  set_up_figure, latex_table, do_table


def subopt_distribution(all_ds, subopts):
  xid = 'strucpersec'
  frames = {}
  for subopt in subopts:
    ds_frames = [ds[subopt][ds[subopt][xid] != float('inf')] for ds in all_ds]
    frames[subopt] = pd.concat(ds_frames)

  do_table(frames, ['strucpersec'], True)
  f, axes = get_subplot_grid(len(frames))
  for i, frame_id in enumerate(sorted(frames.keys())):
    sns.distplot(frames[frame_id][xid], kde=False, bins=50, ax=axes[i], label=frame_id)

  f.tight_layout()
  for ax in axes:
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
  set_up_figure(f, names=(colmap[xid], None), legend=True)
  f.suptitle('Number of structures per second distribution', y=1.00)
  return f


def subopt_perf_results(ds):
  fmap_by_len = {name: frame.groupby('length') for name, frame in ds.fmap.items()}

  # savefig_local(
  #   ds.name, 'strucpersec',
  #   do_quantity_plot(fmap_by_len, 'length', 'strucpersec'))

  # savefig_local(
  #   ds.name, 'real',
  #   do_quantity_plot(fmap_by_len, 'length', 'real'))
  # savefig_local(
  #   ds.name, 'numstruc',
  #   do_quantity_plot(fmap_by_len, 'length', 'numstruc'))
  # savefig_local(
  #   ds.name, 'maxrss',
  #   do_quantity_plot(fmap_by_len, 'length', 'maxrss'))

  # savefig_local(
  #   ds.name, 'usersys',
  #   do_quantity_plot(fmap_by_len, 'length', ['real', 'usersys']))

  # f1, f2 = do_quantity_log_plot(fmap_by_len, 'length', 'real')
  # savefig_local(ds.name, 'real_loglog_scatter', f1)
  # savefig_local(ds.name, 'real_loglog_bestfit', f2)
  #
  # f1, f2 = do_quantity_log_plot(fmap_by_len, 'length', 'maxrss', logx=False)
  # savefig_local(ds.name, 'maxrss_logy_scatter', f1)
  # savefig_local(ds.name, 'maxrss_logy_bestfit', f2)
  #
  # f1, f2 = do_quantity_log_plot(fmap_by_len, 'length', 'maxrss')
  # savefig_local(ds.name, 'maxrss_loglog_scatter', f1)
  # savefig_local(ds.name, 'maxrss_loglog_bestfit', f2)

