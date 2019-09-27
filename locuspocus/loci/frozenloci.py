import uuid
import gzip

from functools import wraps

from locuspocus.locusdb import FrozenLocusDB
from locuspocus.locus import ThawedLocus

from .locimixin import LociMixin


class FrozenLoci(LociMixin,FrozenLocusDB):

    def __init__(self,name,basedir=None):
        super().__init__(name,basedir=basedir)

    def _get_locus_by_UUID(self,UUID):
        return ThawedLocus(UUID,self._db)

    def add_locus(
        self,
        chromosome: str,
        start: int,
        end: int,

        source: str = 'locuspocus',
        feature_type: str = 'locus',
        strand: str = '+',
        frame: int = None,
        name: str = None,

        # Extra locus stuff
        attrs: dict = None,
        parent = None,
        children = None
    ):
        # this starts a transaction
        with self.m80.db.bulk_transaction() as cur:
            # get the fresh UUID
            locus = ThawedLocus(uuid.uuid4().hex,self._db)
            # insert the core feature data
            cur.execute(
                '''
                INSERT INTO loci 
                    (UUID,chromosome,start,end,source,feature_type,strand,frame,name)
                    VALUES (?,?,?,?,?,?,?,?,?)
                ''',(
                # chrom/start/end are required, so cast them
                (locus._UUID,str(chromosome),int(start),int(end),
                 source,feature_type,strand,frame,name))
            )

            # Add the key val pairs
            if attrs is not None:
                for key,val in attrs.items():
                    locus[key] = val          

            # Handle Parent Child Relationships
            locus.parent = parent
            locus.children = children

            return locus
    
    def import_gff(
        self, 
        filename: str, 
        feature_type="*", 
        ID_attr="ID", 
        parent_attr='Parent',
        attr_split="="
    ) -> None:
        '''
            Imports Loci from a gff (General Feature Format) file.
            See more about the format here:
            http://www.ensembl.org/info/website/upload/gff.html

            Parameters
            ----------

            filename : str
                The path to the GFF file.
            ID_attr : str (default: ID)
                The key in the attribute column which designates the ID or
                name of the feature.
            parent_attr : str (default: Parent)
                The key in the attribute column which designates the Parent of
                the Locus
            attr_split : str (default: '=')
                The delimiter for keys and values in the attribute column
        '''
        if filename.endswith(".gz"):
            IN = gzip.open(filename, "rt")
        else:
            IN = open(filename, "r")
        loci = []
        name_index = {}
        total_loci = 0
        for i,line in enumerate(IN):
            total_loci += 1
            # skip comment lines
            if line.startswith("#"):
                continue
            # Get the main information
            (
                chromosome,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                attributes,
            ) = line.strip().split("\t")
            # Cast data into appropriate types
            start = int(start)
            end = int(end)
            strand = None if strand == '.' else strand
            frame = None if frame == '.' else int(frame)
            # Get the attributes
            attributes = dict(
                [
                    (field.strip().split(attr_split))
                    for field in attributes.strip(";").split(";")
                ]
            )
            # Store the score in the attrs if it exists
            if score != '.':
                attributes['score'] = float(score)
            # Parse out the Name (Identifier)
            if ID_attr in attributes:
                name = attributes[ID_attr]
                del attributes[ID_attr]
            else:
                name = None
            # Parse out the parent info
            if parent_attr in attributes:
                parent = self._get_locus_by_UUID(
                    attributes[parent_attr]
                )
                del attributes[parent_attr]
            else:
                parent = None

            # Commit the locus
            locus = self.add_locus(
                chromosome,
                start,
                end,
                source,
                feature_type,
                strand,
                frame,
                name,
                attrs=attributes
            ) 

            if parent is not None:
                locus.parent = parent
        IN.close()
        return total_loci
