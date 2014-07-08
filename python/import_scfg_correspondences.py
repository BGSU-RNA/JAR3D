"""

A program for importing correspondences between motif groups and SCFG models
based on the motif groups. Has to be run every time a new motif group is
released.

Usage:
python import_scfg_correspondences motif_type release_id model_type

Example:
python import_scfg_correspondences IL 0.6 bp_models

NB! relies on a config file called import_scfg.cfg
Template for the config file:
[database]
user:
password:
database:
host:

[paths]
lib: /Users/anton/Dropbox/BGSURNA/Motifs/lib

"""

import os
import re
import sys
import ConfigParser

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.dialects.mysql import LONGTEXT, VARCHAR

import sqlalchemy.exc

# stores database connection info + paths
CONFIG_FILE = 'import_scfg.cfg'

# get command line arguments
motif_type    = sys.argv[1] #'IL'
motif_release = sys.argv[2] #'0.6'
model_type    = sys.argv[3] #'bp_models'
group_set     = motif_type + motif_release

# import configs
config = ConfigParser.RawConfigParser()
config.read(CONFIG_FILE)
user     = config.get('database', 'user')
password = config.get('database', 'password')
host     = config.get('database', 'host')
database = config.get('database', 'database')
lib      = config.get('paths', 'lib') # must be an OS-specific full path

# database setup
engine  = create_engine('mysql://'+user+':'+password+'@'+host+'/'+database)
Session = sessionmaker(bind=engine)
session = Session()

Base = declarative_base()

# describe the jar3d_correspondences table
class SCFGCorrespondence(Base):
    """
    """
    __tablename__ = 'jar3d_correspondences'

    id         = Column(Integer, primary_key=True) # django primary key
    motif_id   = Column(String(11)) # IL_12345.XX
    group_set  = Column(String(6))  # IL0.6
    model_type = Column(String(20)) # bp_models
    correspondences = Column(Text)  # correspondences

Base.metadata.create_all(engine)

# import files
data_folder = os.path.join(lib,
                           motif_type,
                           motif_release,
                           model_type)

files = os.listdir(data_folder)

for file in files:
    if re.search('_correspondences.txt', file):
        print file
        f = open( os.path.join(data_folder, file) )

        m = re.match("(.+)_correspondences.txt", file)
        motif_id = m.group(1)
        correspondences = f.read()

        session.add( SCFGCorrespondence(motif_id = motif_id,
                                        group_set = group_set,
                                        model_type = model_type,
                                        correspondences = correspondences) )
session.commit()

print 'Done'