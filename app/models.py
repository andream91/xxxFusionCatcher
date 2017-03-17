from neomodel import (StructuredNode, StructuredRel, StringProperty, IntegerProperty, ArrayProperty)
from neomodel.relationship_manager import RelationshipTo, RelationshipFrom
from neomodel.properties import FloatProperty

#EDGES
class AT_CHROMOSOME(StructuredRel):
    fusion_partner = StringProperty()

class IN_COUPLE(StructuredRel):
    position = IntegerProperty()
    
    
class WITH_OTHER_TRANSCRIPT(StructuredRel):
    position = IntegerProperty()   
    
class WITH_VIRUSES(StructuredRel):
    count_of_mapping_reads = IntegerProperty()

#NODES
class CellLine(StructuredNode):
    cell_line = StringProperty()
    #
    happen = RelationshipTo('Fusion',"HAPPEN")
    with_viruses = RelationshipTo('Virus',"WITH_VIRUSES",model=WITH_VIRUSES)
    
class Chromosome(StructuredNode):
    chromosome = StringProperty()
    #
    of_gene = RelationshipTo('Gene',"OF_GENE")
    #
    fromFusiontoChromosome = RelationshipFrom('Fusion', "AT_CHROMOSOME", model=AT_CHROMOSOME)
    
class Couple(StructuredNode):
    couple = IntegerProperty() 
    #
    with_other_transcript = RelationshipTo('Transcript',"WITH_OTHER_TRANSCRIPT", model=WITH_OTHER_TRANSCRIPT)
    with_protein = RelationshipTo('Protein',"WITH_PROTEIN")
    #
    fromTranscriptToCouple = RelationshipFrom('Transcript',"IN_COUPLE", model=IN_COUPLE)
    fromFusionToCouple = RelationshipFrom('FusionCatcher',"WITH_TRANS_COUPLE")
    
class EricScript(StructuredNode):
    ericscript_id = IntegerProperty()
    breakpoint_1 = StringProperty()
    strand_1 = StringProperty()
    breakpoint_2 = StringProperty()
    strand_2 = StringProperty()
    crossing_reads = IntegerProperty()
    spanning_reads = IntegerProperty()
    mean_intersize = FloatProperty()
    homology = StringProperty()
    fusion_type = StringProperty()
    junction_sequence = StringProperty()
    gene_expr_1 = FloatProperty()
    gene_expr_2 = FloatProperty()
    gene_expr_fused = FloatProperty()
    es = FloatProperty()
    gjs = StringProperty()
    us = FloatProperty()
    eric_score = FloatProperty()

    #
    fromFusionToEricScript = RelationshipFrom('Fusion', "WITH_ERIC_SCRIPT")
    
class Exon(StructuredNode):
    exon = StringProperty()
    #
    in_gene = RelationshipTo('Gene',"IN_GENE")
    #
    fromFusionToExon = RelationshipFrom('FusionCatcher', "AT_EXON")
    
class Fusion(StructuredNode):
    fusion_id = IntegerProperty()
    #
    with_fc_script = RelationshipTo('FusionCatcher',"WITH_FC_SCRIPT")
    with_eric_script =  RelationshipTo('EricScript',"WITH_ERIC_SCRIPT")
    with_tophat_script = RelationshipTo('Tophat',"WITH_TOPHAT_SCRIPT")
    with_gene = RelationshipTo('Gene',"WITH")
    at_chromosome = RelationshipTo('Chromosome',"AT_CHROMOSOME")
    #
    fromCellLineToFusion = RelationshipFrom('CellLine',"HAPPEN")
    fromGeneToFusion = RelationshipFrom('Gene',"HAD")
    
    
class FusionCatcher(StructuredNode):
    fusion_id = IntegerProperty()
    description = ArrayProperty()
    common_mapping_reads = IntegerProperty()
    spanning_pairs = IntegerProperty()
    spanning_unique_reads = IntegerProperty()
    longest_anchor_found = IntegerProperty()
    fusion_finding_method = StringProperty()
    fusion_sequence = StringProperty()
    fusion_point_1 = IntegerProperty()
    fusion_point_2 = IntegerProperty()
    strand_1 = StringProperty()
    strand_2 = StringProperty()
    predicted_effect_1 = StringProperty()
    predicted_effect_2 = StringProperty()
    #
    at_exon = RelationshipTo('Exon', "AT_EXON")
    with_trans_couple = RelationshipTo('Couple', "WITH_TRANS_COUPLE")
    with_gene = RelationshipTo('Gene',"WITH")
    #
    #fromCellLineToFusion = RelationshipFrom('CellLine',"HAPPEN")
    fromFusionToFusionCatcher = RelationshipFrom('Fusion', "WITH_FC_SCRIPT")

class Gene(StructuredNode):
    ensid = StringProperty()
    symbol = StringProperty()
    description = StringProperty()
    #
    had = RelationshipTo('Fusion', "HAD")
    #
    #fromFusion = RelationshipFrom('Fusion',"WITH", model=WITH)
    fromExonToGene = RelationshipFrom('Exon',"IN_GENE")
    fromChromosomeToGene = RelationshipFrom('Chromosome',"OF_GENE")
    fromFusionToGene = RelationshipFrom('Fusion','WITH')
    
class Protein(StructuredNode):
    protein = StringProperty()
    #
    fromCoupleToProtein = RelationshipFrom('Couple',"WITH_PROTEIN")
    
class Tophat(StructuredNode):
    tophat_id = IntegerProperty()
    left_coord = StringProperty()
    right_coord = StringProperty()
    spanning_reads = StringProperty()
    spanning_mate_pairs = StringProperty()
    spanning_mate_pairs_end = StringProperty()
    #
    fromFusionToTophat = RelationshipFrom('Fusion','WTIH_TOPHAT_SCRIPT')
    
class Transcript(StructuredNode):
    transcript = StringProperty()
    #
    in_couple = RelationshipTo('Couple',"IN_COUPLE", model=IN_COUPLE)
    #
    fromCoupleToTranscript = RelationshipFrom('Couple',"WITH_OTHER_TRANSCRIPT", model=WITH_OTHER_TRANSCRIPT)

    
class Virus(StructuredNode):
    name = StringProperty()
    gi = StringProperty()
    ref = StringProperty()
    #
    fromCellLineToVirus = RelationshipFrom('CellLine',"WITH_VIRUSES", model=WITH_VIRUSES)