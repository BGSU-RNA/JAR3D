import sys
import os
sys.path.append('/Users/api/apps/jar3d_dev/app/JAR3Doutput')  #set this to jar3d-website folder 
from myproject import settings
from django.core.management import setup_environ
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "JAR3Dresults.settings")


from JAR3Dresults.models import Query_info
from JAR3Dresults.models import Query_sequences
from JAR3Dresults.models import Query_loop_positions
from JAR3Dresults.models import Results_by_loop

from JAR3Dresults.views import make_query_sequences
from JAR3Dresults.views import make_query_indices
from JAR3Dresults.views import make_query_info
from JAR3Dresults.views import format_extracted_loops

from JAR3Dresults.settings import DATABASES

from Bio import SeqIO
import subprocess

def main(arg):
	query_id = str(uuid.uuid4())
	infile = open(arg.infile)
	lines = infile.readlines()
	ss = lines[0]
	data = lines[2::2]
	fasta = lines[1::2]
	parsed_input = infile.read()
	infile.close()
	for index, item in lines:
		if index
	query_type = 'isFastaMultipleSequencesSS'
	try:
        loops, indices = self.isfolded_extract_loops(ss, data)
        except fold.FoldingTimeOutError:
            return self.respond("Folding timed out")
        except fold.FoldingFailedError:
            return self.respond("Folding failed")
        except:
            return self.respond("Couldn't extract loops")

    query_info = JAR3Dresults.views.make_query_info(query_id, query_type, parsed_input)
    query_sequences = JAR3Dresults.views.make_query_sequences(loops, fasta, query_id)
    query_positions = JAR3Dresults.views.make_query_indices(indices, query_id)

    # don't proceed unless there are internal loops
    if not query_sequences:
        sys.exit("No internal loops found in the input")

    # todo: if all loops have status = -1, then set query_info.status to 1
    # persist the entries in the database starting with sequences
    # query_sequences, query_positions, mins = zip(*sort_sequences(query_positions, query_sequences))
    try:
        for seq in query_sequences:
            seq.save()
    except:
        sys.exit("Couldn't save query_sequences")
    try:
        for ind in query_positions:
            ind.save()
    except:
        sys.exit("Couldn't save query_positions")
    try:
        query_info.save()
    except:
        sys.exit("Couldn't save query_info")
    # everything went well, pass off to jar3d jar to run query
    IL_path = "/Users/api/Models/IL/1.13/lib/all.txt"
    HL_path = "/Users/api/Models/HL/1.13/lib/all.txt"
    db_name = JAR3Dresults.settings.DATABASES['default']['NAME']
    db_usr_name = JAR3Dresults.settings.DATABASES['default']['USER']
    db_pw = JAR3Dresults.settings.DATABASES['default']['PASSWORD']
    jar_path = "/Users/api/apps/jar3d_dev/app/queue/webJAR3D_server.jar"
    cmd = ['java', '-jar', jar_path, IL_path, HL_path, query_id, db_usr_name, db_pw, db_name]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	for line in p.stdout:
		print line
	p.wait()

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('infile', help='Path to Fasta input file')
    arg = p.parse_args()
    main(arg)