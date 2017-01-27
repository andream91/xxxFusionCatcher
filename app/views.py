from app.models import CellLine, Chromosome, Gene, Fusion, Exon, Transcript
import json, csv
from neomodel import db, match

# Create your views here.
from django.http import HttpResponse


def search_for_cell_line(request,c_line):
    response = {}
    header = ["Cell line",
        "Gene pair symbols",
        "Gene pair EnsIDs",
        "Exon pair",
        "Chromosome : fusion point : strand",
        "Description",
        "Counts of common mapping reads",
        "Spanning pairs",
        "Spanning unique reads",
        "Longest anchor found",
        "Fusion finding method",
        "Fusion sequence",
        "Predicted effect",
        "Predicted fused transcripts",
        "Predicted fused proteins"]
    
    # recupero fusioni nella linea cellulare
    fusions = []
    #se per tutte le linee cellulari mostra tutte le fusioni
    if c_line == "ALL":
        for c_l in CellLine.nodes.all():
            for fusion in c_l.happen:
                fusions.append(fusion)
                #print(fusion)
    else:
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)
    
    rows = build_rows(fusions,header)
    #print(rows)
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_chromosome(request,c_line,chromos,start_point,end_point):
    response = {}
    header = ["Cell line",
        "Gene pair symbols",
        "Gene pair EnsIDs",
        "Exon pair",
        "Chromosome : fusion point : strand",
        "Description",
        "Counts of common mapping reads",
        "Spanning pairs",
        "Spanning unique reads",
        "Longest anchor found",
        "Fusion finding method",
        "Fusion sequence",
        "Predicted effect",
        "Predicted fused transcripts",
        "Predicted fused proteins"]
    
    # recupero fusioni nella linea cellulare
    fusions = []
    #se ho specificato solo il cromosoma, cerco tutte le fusioni in tutte le linee cellulari che coinvolgono il cromosoma
    if c_line == "ALL" and chromos != "":
        c = Chromosome.nodes.get(chromosome = chromos)
        for fusion in c.fromFusiontoChromosome:
            if fusion.fusion_point_1 >= start_point and fusion.fusion_point_2 <= end_point:
                fusions.append(fusion)
    #se ho specificato solo la linea cellulare, cerco tutte le fusioni che avvengono in uqella linea cellulare nel determinato intervallo
    elif c_line != "" and chromos == "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.fusion_point_1 >= start_point and fusion.fusion_point_2 <= end_point:
                fusions.append(fusion)
    #se ho sia linea cellulare che cromosoma specificati, cerco tutte le fusioni nella linea cellulare che coinvolgono il cromosoma
    elif c_line != "" and chromos != "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.fusion_point_1 >= int(start_point) and fusion.fusion_point_2 <= int(end_point):
                if fusion.at_chromosome.filter(chromosome__exact=chromos):
                    fusions.append(fusion)
    
    rows = build_rows(fusions,header)
    print(rows)
    
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_gene(request,c_line,gene_one,gene_two):
    response = {}
    header = ["Cell line",
        "Gene pair symbols",
        "Gene pair EnsIDs",
        "Exon pair",
        "Chromosome : fusion point : strand",
        "Description",
        "Counts of common mapping reads",
        "Spanning pairs",
        "Spanning unique reads",
        "Longest anchor found",
        "Fusion finding method",
        "Fusion sequence",
        "Predicted effect",
        "Predicted fused transcripts",
        "Predicted fused proteins"]
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #se ho specificato il gene ma non la linea cellulare, mostro tutte le fusioni che coinvolgono il gene
    if c_line == "ALL" and gene_one != "" and gene_two == "":
        if "ENSG" in gene_one:
            g = Gene.nodes.get(gene_id = gene_one)
            if g.had:
                for fusion in g.had:
                    fusions.append(fusion)
            if g.fromFusionToGene:
                for fusion in g.fromFusionToGene:
                    fusions.append(fusion)
        else:
            g = Gene.nodes.get(symbol = gene_one)
            if g.had:
                for fusion in g.had:
                    fusions.append(fusion)
            if g.fromFusionToGene:
                for fusion in g.fromFusionToGene:
                    fusions.append(fusion)
    elif c_line == "ALL" and gene_one != "" and gene_two != "":
        if "ENSG" in gene_one and "ENSG" in gene_two:
            g = Gene.nodes.get(gene_id = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].gene_id == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(gene_id = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].gene_id == gene_one:
                        fusions.append(fusion)
        elif "ENSG" in gene_one and "ENSG" not in gene_two:
            g = Gene.nodes.get(gene_id = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(symbol = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].gene_id == gene_one:
                        fusions.append(fusion)
        elif "ENSG" not in gene_one and "ENSG" in gene_two:
            g = Gene.nodes.get(symbol = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].gene_id == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(gene_id = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_one:
                        fusions.append(fusion)
        else:
            g = Gene.nodes.get(symbol = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_two:
                        fusions.append(fusion)
                        
            g = Gene.nodes.get(symbol = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_one:
                        fusions.append(fusion)    
    #se ho specificato la linea cellulare ma non il gene, mostro tutte le fusioni per quella linea cellulare (ANALOGO A SEARCH FOR CELL_LINE)
    elif c_line != "ALL" and gene_one == "" and gene_two == "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)     
    #se ho specificato sia la linea cellulare che il gene, cerco tutte le fusioni che coinvolgono quel gene nella linea cellulare
    elif c_line != "ALL" and gene_one != "" and gene_two == "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.fromGeneToFusion.filter(symbol__exact=gene_one) or fusion.fromGeneToFusion.filter(gene_id__exact=gene_one) or fusion.with_gene.filter(symbol__exact=gene_one) or fusion.with_gene.filter(gene_id__exact=gene_one):
                fusions.append(fusion)
    elif c_line != "ALL" and gene_one != "" and gene_two != "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if (fusion.fromGeneToFusion.filter(symbol__exact=gene_one) and fusion.with_gene.filter(symbol__exact=gene_two)) or (fusion.fromGeneToFusion.filter(symbol__exact=gene_two) and fusion.with_gene.filter(symbol__exact=gene_one)) or (fusion.fromGeneToFusion.filter(gene_id__exact=gene_one) and fusion.with_gene.filter(gene_id__exact=gene_two)) or (fusion.fromGeneToFusion.filter(gene_id__exact=gene_two) and fusion.with_gene.filter(gene_id__exact=gene_one)):
                fusions.append(fusion)
            
    rows = build_rows(fusions,header)
    print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))



def search_for_exon(request,c_line,exon_one,exon_two):
    response = {}
    header = ["Cell line",
        "Gene pair symbols",
        "Gene pair EnsIDs",
        "Exon pair",
        "Chromosome : fusion point : strand",
        "Description",
        "Counts of common mapping reads",
        "Spanning pairs",
        "Spanning unique reads",
        "Longest anchor found",
        "Fusion finding method",
        "Fusion sequence",
        "Predicted effect",
        "Predicted fused transcripts",
        "Predicted fused proteins"]
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #cell_line all, primo esone si, secondo esone no, tutte le fusioni che coinvolgono questo esone
    if c_line == "ALL" and exon_one!="" and exon_two == "":
        e = Exon.nodes.get(exon = exon_one)
        for fusion in e.fromFusionToExon:
            fusions.append(fusion)
    #cell_line all, primo esone si, secondo esone si (011), tutte le fusioni che convolgono la coppia di esoni
    elif c_line == "ALL" and exon_one!="" and exon_two!="":
        e = Exon.nodes.get(exon = exon_one)
        for fusion in e.fromFusionToExon:
            if fusion.at_exon.filter(exon__exact=exon_two):
                fusions.append(fusion)
    #cell_line si, primo esone no, secondo esone no, ANALOGO A SEARCH FOR CELL_LINE
    elif c_line!="ALL" and exon_one=="" and exon_two=="":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)     
    #cell_line si, primo esone si, secondo esone no 
    elif c_line != "ALL" and exon_one != "" and exon_two == "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.at_exon.filter(exon__exact=exon_one):
                fusions.append(fusion)
    #cell_line si, primo esone si, secondo esone si 
    elif c_line != "ALL" and exon_one!="" and exon_two != "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.at_exon.filter(exon__exact=exon_one) and fusion.at_exon.filter(exon__exact=exon_two):
                    fusions.append(fusion)

    rows = build_rows(fusions,header)
    print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_transcript(request,c_line,transcript_one,transcript_two):
    response = {}
    header = ["Cell line",
        "Gene pair symbols",
        "Gene pair EnsIDs",
        "Exon pair",
        "Chromosome : fusion point : strand",
        "Description",
        "Counts of common mapping reads",
        "Spanning pairs",
        "Spanning unique reads",
        "Longest anchor found",
        "Fusion finding method",
        "Fusion sequence",
        "Predicted effect",
        "Predicted fused transcripts",
        "Predicted fused proteins"]
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #linea cellulare non specificata, un trascritto specificato
    if c_line=="ALL" and transcript_one!="" and transcript_two=="":
        for couple in Transcript.nodes.get(transcript = transcript_one).fromCoupleToTranscript:
            for fusion in couple.fromFusionToCouple:
                fusions.append(fusion)
    #linea cellulare non specificata, coppia di trascritti specificata
    elif c_line=="ALL" and transcript_one!="" and transcript_two!="":
        for couple in Transcript.nodes.get(transcript = transcript_one).fromCoupleToTranscript:
            if couple.with_other_transcript.filter(transcript__exact=transcript_two) or couple.fromTranscriptToCouple.filter(transcript__exact=transcript_two):
                for fusion in couple.fromFusionToCouple:
                    fusions.append(fusion)
        for couple in Transcript.nodes.get(transcript = transcript_two).fromCoupleToTranscript:
            if couple.with_other_transcript.filter(transcript__exact=transcript_one) or couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one):
                for fusion in couple.fromFusionToCouple:
                    fusions.append(fusion)
    #linea cellulare specificata, un trascritto specificato:
    elif c_line!="ALL" and transcript_one!="" and transcript_two=="trolo":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            for couple in fusion.with_trans_couple:
                if couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one) or couple.with_other_transcript.filter(transcript__exact=transcript_one):
                    fusions.append(fusion)
    #linea cellulare specificata, coppia di trascritti specificata
    elif c_line!="ALL" and transcript_one!="" and transcript_two!="":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            for couple in fusion.with_trans_couple:
                if (couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one) and couple.with_other_transcript.filter(transcript__exact=transcript_two)) or (couple.fromTranscriptToCouple.filter(transcript__exact=transcript_two) and couple.with_other_transcript.filter(transcript__exact=transcript_one)):
                    fusions.append(fusion)
        

    rows = build_rows(fusions,header)
    print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_fusion_information(request,c_line,algorithm,fusion_description,predicted_effect1,predicted_effect2):
    response = {}
    header = ["Cell line",
        "Gene pair symbols",
        "Gene pair EnsIDs",
        "Exon pair",
        "Chromosome : fusion point : strand",
        "Description",
        "Counts of common mapping reads",
        "Spanning pairs",
        "Spanning unique reads",
        "Longest anchor found",
        "Fusion finding method",
        "Fusion sequence",
        "Predicted effect",
        "Predicted fused transcripts",
        "Predicted fused proteins"]
    
    # recupero fusioni nella linea cellulare
    fusions = []
    for fusion in CellLine.nodes.get(cell_line = c_line).happen:
        predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
        predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
        if (algorithm in fusion.fusion_finding_method) and (fusion_description in fusion.description) and (predicted_effect1 == predicted_effect_1) and (predicted_effect2 == predicted_effect_2):
            fusions.append(fusion)

    rows = build_rows(fusions,header)
    print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def build_rows(fusions, header):
    rows = []
    # ora che ho solo le fusioni interessate recupero le informazioni e mi costruisco la riga
    for myfusion in fusions:
        # recupero cell line
        cellLine = myfusion.fromFusionToCellLine.all()[0].cell_line
        
        #recupero dati dai geni
        gene1 = myfusion.fromGeneToFusion.all()[0]
        strand_1 = myfusion.fromGeneToFusion.relationship(gene1).strand
        predicted_effect_1 = myfusion.fromGeneToFusion.relationship(gene1).predicted_effect
        gene2 = myfusion.with_gene.all()[0]
        strand_2 = myfusion.with_gene.relationship(gene2).strand
        predicted_effect_2 = myfusion.with_gene.relationship(gene2).predicted_effect

        #recupero cromosomi 
        chromosome1 = myfusion.at_chromosome.match(fusion_partner__exact="5'end").all()[0]
        chromosome2 = myfusion.at_chromosome.match(fusion_partner__exact="3'end").all()[0]
        fusion_point_1 = myfusion.fusion_point_1
        fusion_point_2 = myfusion.fusion_point_2

        #recupero esoni        
        exon1 = []
        exon2 = []

        for exon in myfusion.at_exon:
            if exon.in_gene.filter(symbol__exact=gene1.symbol):
                exon1 = exon
            if exon.in_gene.filter(symbol__exact=gene2.symbol):
                exon2 = exon
                
        #recupero trascritti e proteine
        transcript_couples = []
        proteins = []
        for couple in myfusion.with_trans_couple:
            transcript1 = couple.fromTranscriptToCouple.all()[0]
            transcript1_position = couple.fromTranscriptToCouple.relationship(transcript1).position
            transcript2 = couple.with_other_transcript.all()[0]
            transcript2_position = couple.with_other_transcript.relationship(transcript2).position
            
            transcript_couples.append(transcript1.transcript+":"+str(transcript1_position)+" - "+transcript2.transcript+":"+str(transcript2_position))
            if couple.with_protein.all():
                proteins.append(couple.with_protein.all()[0].protein)
            
        #costruisco la riga
        row = []
        row.append(cellLine)
        row.append(gene1.symbol+" - "+gene2.symbol)
        row.append(gene1.gene_id+" - "+gene2.gene_id)
        if exon1 or exon2:
            row.append(exon1.exon+" - "+exon2.exon)
        else:
            row.append("No exons")
        row.append([chromosome1.chromosome+":"+str(fusion_point_1)+":"+strand_1, chromosome2.chromosome+":"+str(fusion_point_2)+":"+strand_2])
        row.append(myfusion.description)
        row.append(myfusion.common_mapping_reads)
        row.append(myfusion.spanning_pairs)
        row.append(myfusion.spanning_unique_reads)
        row.append(myfusion.longest_anchor_found)
        row.append(myfusion.fusion_finding_method)
        row.append(myfusion.fusion_sequence)
        if(predicted_effect_1!=predicted_effect_2):
            row.append(predicted_effect_1+"/"+predicted_effect_2)
        else:
            row.append(predicted_effect_1)
        row.append(transcript_couples)
        row.append(proteins)   
        #print(row)
        rows.append(row)
    return rows

def search_viruses(request,c_line,vir):
    response = {}
    rows = []
    header = ["Cell line",
        "Virus/bacteria/phage name",
        "GI",
        "NC"]
    
    #recupero virus nella linea cellulare
    #se ho specificato solo il virus, cerca per tutte le linee cellulari
    if c_line == "ALL" and vir != "":
        for c_l in CellLine.nodes.all():
            for virus in c_l.with_viruses:
                if vir in virus.name or vir == virus.gi or vir == virus.ref:
                    rows.append([c_l.cell_line, virus.name, virus.gi, virus.ref])
    #se ho specificato solo la linea cellulare, cerca tutti i virus per la linea cellulare
    elif c_line !="" and vir == "trolo":
        for virus in CellLine.nodes.get(cell_line = c_line).with_viruses:
            rows.append([c_line, virus.name, virus.gi, virus.ref])
    #se ho specificato sia virus che linea cellulare, cerca quel virus per quella linea cellulare
    elif c_line != "" and vir != "":
        for virus in CellLine.nodes.get(cell_line = c_line).with_viruses:
            if vir in virus.name or vir == virus.gi or vir == virus.ref:
                rows.append([c_line, virus.name, virus.gi, virus.ref])
                
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def generate_statistics(request):
    #chromosome-fusion
    chromosome_fusion_f =  open('chromosome_fusion.csv','w')
    chromosome_fusion_w = csv.writer(chromosome_fusion_f)
    chromosome_fusion_w.writerow(["Chromosome","Fusion"])
    for chromosome in Chromosome.nodes.all():
        #print(len(chromosome.fromFusiontoChromosome))
        chromosome_fusion_w.writerow([chromosome.chromosome,len(chromosome.fromFusiontoChromosome)])
    chromosome_fusion_f.close() 
    
    #cell line-fusion
    cell_line_fusion_f =  open('cell_line_fusion.csv','w')
    cell_line_fusion_w = csv.writer(cell_line_fusion_f)
    cell_line_fusion_w.writerow(["Cell Line","Fusion"])
    for cell_line in CellLine.nodes.all():
        cell_line_fusion_w.writerow([cell_line.cell_line,len(cell_line.happen)])
    cell_line_fusion_f.close()

    #cell line-gene
    cell_line_gene_f =  open('cell_line_gene.csv','w')
    cell_line_gene_w = csv.writer(cell_line_gene_f)
    cell_line_gene_w.writerow(["Cell Line","Gene"])
    #for cell_line in CellLine.nodes.all():
    #   query = "match (c:CellLine{cell_line:'"+cell_line.cell_line+"'})-[:HAPPEN]->(f:Fusion) with c, f match (g1:Gene)-[:HAD]->(f) WITH c,f,collect(DISTINCT g1) AS set1 match (f)-[:WITH]->(g2:Gene) with c,f,set1,collect(DISTINCT g2) AS set2 with set1+set2 as both unwind both as res return count(distinct res)"
    #   results = db.cypher_query(query)
    #   cell_line_gene_w.writerow([cell_line.cell_line,db.cypher_query(query)[0][0][0]])
    #for cell_line in CellLine.nodes.all():
    #    genes = []
    #    for fusion in cell_line.happen:
    #        if fusion.fromGeneToFusion.all()[0].symbol not in genes:
    #                 genes.append(fusion.fromGeneToFusion.all()[0].symbol)
    #             if fusion.with_gene.all()[0].symbol not in genes:
    #                 genes.append(fusion.with_gene.all()[0].symbol)
    #         cell_line_gene_w.writerow([cell_line.cell_line,len(genes)])
    #         print("done "+cell_line.cell_line)
    for cell_line in CellLine.nodes.all():
        definition = dict(node_class=Gene, direction=match.OUTGOING, relation_type='*', model=None)
        relations_traversal = match.Traversal(cell_line, Gene.__label__, definition)
        genes = relations_traversal.all()
        print(genes)
        
    cell_line_gene_f.close()
    
    return HttpResponse()
        
    
    
    
    
    
    
    