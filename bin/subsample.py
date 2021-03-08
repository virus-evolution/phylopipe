class Node:
    def __init__(self, i, mutations=[], samples=[], outnodes=[], parent_id=None):
        self.node_id = i
        self.mutations = mutations.copy()
        self.samples = samples.copy()
        self.outnodes = outnodes.copy()
        self.parent_id = parent_id
        self.fuzzy = False
        
    def intersection(self,mutation_list):
        return [x for x in self.mutations if x in mutation_list]
        
    def relationship(self,intersection,fuzzy=0):
        if len(intersection) == 0:
            return "unrelated"
        elif len(intersection) == len(self.mutations):
            return "descendant"
        elif len(intersection) > 0 and len(intersection) >= len(self.mutations) - fuzzy:
            return "fuzzy_descendant"
        else:
            return "partial"
        
    def update_on_split(self,mutations,parent_id):
        self.mutations = mutations.copy()
        self.parent_id = parent_id
        #print("updating node:",self.node_id,self.mutations, self.samples, self.outnodes, self.parent_id)
        
    def update_on_prune(self,mutations,samples,outnode):
        self.mutations.extend([x for x in mutations if x not in self.mutations])
        self.samples.extend(samples)
        self.outnodes.remove(outnode)
        self.fuzzy = True
        #print("updating node:",self.node_id, self.mutations, self.samples, self.outnodes, self.parent_id)
        
    def write(self, verbose):
        if verbose:
            print("Node:", self.node_id, "mutations:", self.mutations, "samples:", self.samples, "outnodes:", self.outnodes, "parent:", self.parent_id)
        else:
            print("Node:", self.node_id, "samples:", self.samples, "outnodes:", self.outnodes, "parent:", self.parent_id)
    
    def to_tsv(self):
        return "\t".join(str(self.node_id), str(self.mutations), str(self.samples), str(self.outnodes), str(self.parent_id), str(self.fuzzy))
        

class Graph:
    def __init__(self):
        self.next_i = 0
        self.root = 0
        self.nodes = {}
        self.ordered_nodes = []
        self.add_node()
        
    def add_node(self,mutations=[], samples=[], outnodes=[], parent_id=None):
        #print("adding node:",self.next_i,mutations, samples, outnodes, parent_id)
        self.nodes[self.next_i] = Node(self.next_i,mutations, samples, outnodes, parent_id)
        if parent_id != None:
            self.nodes[parent_id].outnodes.append(self.next_i)
        self.ordered_nodes.append(self.next_i)
        self.next_i += 1
        
        
    def split_node(self, node_id, sample, mutation_list):
        #print("split node", node_id, "for sample", sample)
        left_split = [x for x in self.nodes[node_id].mutations if x not in mutation_list]
        right_split = [x for x in mutation_list if x not in self.nodes[node_id].mutations]
        intersection = [x for x in mutation_list if x in self.nodes[node_id].mutations]
        self.nodes[self.nodes[node_id].parent_id].outnodes.remove(node_id)
        self.add_node(intersection,parent_id=self.nodes[node_id].parent_id, outnodes=[node_id])
        self.ordered_nodes.insert(self.ordered_nodes.index(node_id),self.next_i-1)
        self.ordered_nodes.pop()
        #print(self.ordered_nodes)
        if node_id == self.root:
            self.root = self.next_i-1
        self.nodes[node_id].update_on_split(intersection+left_split,self.next_i-1)
        #self.write()
        if len(right_split) > 0:
            self.add_node(intersection+right_split,[sample],parent_id=self.next_i-1)
        else:
            self.nodes[self.next_i-1].samples.append(sample)
        #self.write()

    
    def find_descendant(self, mutations, fuzzy):
        intermediate_nodes = [self.root]
        final_node = self.root
        final_intersection = 0
        
        def process_intermediate_node():
            nonlocal intermediate_nodes, final_node, final_intersection
            node = self.nodes[intermediate_nodes.pop()]
            
            if len(node.outnodes) == 0:
                if len(node.intersection(mutations)) > final_intersection:
                        final_node = node.node_id
                        final_intersection = len(node.intersection(mutations))
                
            for outnode_id in node.outnodes:
                outnode = self.nodes[outnode_id]
                intersection = outnode.intersection(mutations)
                relationship = outnode.relationship(intersection, fuzzy)
                if relationship == "partial" or (relationship.endswith("descendant") and len(outnode.outnodes) == 0):
                    if len(intersection) > final_intersection:
                        final_node = outnode_id
                        final_intersection = len(intersection)
                elif relationship.endswith("descendant"):
                    intermediate_nodes.append(outnode_id)
                    
        while len(intermediate_nodes) > 0:
            process_intermediate_node()
                    
        return final_node, final_intersection
            
    def add_sample(self, sample, mutations, fuzzy=0):
        #print("Add sample", sample)
        old_num_samples = self.number_sample_classes()
        
        candidate_node, candidate_intersection_size = self.find_descendant(mutations, fuzzy)
        #print(candidate_node, candidate_intersection_size)
        
        if candidate_node == self.root:
            #print("add new node for sample", sample)
            self.add_node(mutations,[sample],parent_id=self.root)
        elif candidate_intersection_size == len(self.nodes[candidate_node].mutations) == len(mutations):
            #print("append sample", sample)
            self.nodes[candidate_node].samples.append(sample)
        elif candidate_intersection_size >= len(self.nodes[candidate_node].mutations) - fuzzy and candidate_intersection_size >= len(mutations) - fuzzy:
            #print("append sample", sample)
            self.nodes[candidate_node].samples.append(sample)
            if not candidate_intersection_size == len(self.nodes[candidate_node].mutations) == len(mutations):
                self.nodes[candidate_node].fuzzy = True
        elif candidate_intersection_size >= len(self.nodes[candidate_node].mutations) - fuzzy:
            #print("add sample with a parent", sample)
            self.add_node(mutations,[sample],parent_id=candidate_node)
            if candidate_intersection_size != len(self.nodes[candidate_node].mutations):
                self.nodes[self.next_i - 1].fuzzy = True
        else:
            #print("split node for sample", sample)
            self.split_node(candidate_node,sample,mutations)
            
        if self.number_sample_classes() > old_num_samples + 1:
            print("Problem adding sample", sample)
            self.debug_sample_classes()            
            
    def number_sample_classes(self):
        count = 0
        for node_id in self.nodes:
            if len(self.nodes[node_id].samples) > 0:
                count += 1
            #print(count, node_id, self.nodes[node_id].samples)
        return count            
        
    def prune(self,distance=None, sample_size=None, internal=False):
        del_keys = []
        for node_id in reversed(self.ordered_nodes):
            node = self.nodes[node_id]
            if node.parent_id is None:
                continue
            elif not internal and len(node.outnodes) > 0:
                continue
            elif sample_size and len(node.samples) > sample_size:
                continue
            elif distance is None:
                del_keys.append(node_id)
            elif len(node.mutations) - len(self.nodes[node.parent_id].mutations) <= distance:
                #print(len(node.mutations),"-",len(self.nodes[node.parent_id].mutations),"<=", distance)
                #print("delete", node_id)
                del_keys.append(node_id)
                
        whitelist = []
        for node_id in del_keys:
            #print("remove", node_id, "and update parent", node.parent_id, "whitelist:", whitelist)
            if node_id in whitelist:
                #print("cancelled - node in whitelist")
                continue
            node = self.nodes[node_id]
            self.nodes[node.parent_id].update_on_prune(node.mutations, node.samples, node_id)
            whitelist.append(node.parent_id)
            del self.nodes[node_id] 
            self.ordered_nodes.remove(node_id)
                       
    def write(self, verbose=False):
        for node_id in self.nodes:
            #print(node_id, "->", self.nodes[node_id].outnodes)
            self.nodes[node_id].write(verbose)
    
    def debug_sample_classes(self):
        count = 0
        for node_id in self.nodes:
            if len(self.nodes[node_id].samples) > 0:
                count += 1
            print(count, node_id, self.nodes[node_id].samples)
            
    def output_sample_classes(self, outfile):
        with open(outfile, "w") as f:
            for node_id in self.nodes:
                for sample in self.nodes[node_id].samples:
                    f.write("%s,%i\n" %(sample, node_id))

    def to_tsv(self, tsv_file):
        with open(tsv_file,"w") as out:
            for node in self.nodes:
                out.write(node.to_tsv())


import csv

in_metadata = "/Users/rmnorris/Work/git/grapevine-global-tree/test/cog_global_2021-02-23_mutations.csv" 
sample_dict = {}
seen = []
with open(in_metadata,"r") as f:
    reader = csv.DictReader(f)
    for r in reader:
        if len(r["variants"]) > 0 and r["variants"] not in seen:
            seen.append(r["variants"])
            sample_dict[r["sequence_name"]] = r["variants"].split("|")
sample_dict = {key: value for key, value in sorted(sample_dict.items(), key=lambda item: item[1])}
print("sample_dict has", len(sample_dict), "keys")

graph = Graph()
fuzzy = 0
count = 0

for sample in sample_dict:
    graph.add_sample(sample,sample_dict[sample],fuzzy)    
    count += 1
    if graph.number_sample_classes() > count:
        graph.debug_sample_classes()
        break
print(graph.number_sample_classes())
graph.to_tsv("/Users/rmnorris/Work/git/grapevine-global-tree/test/cog_global_2021-02-23_mutations.fuzzy0.tsv")
graph.output_sample_classes("/Users/rmnorris/Work/git/grapevine-global-tree/test/cog_global_2021-02-23_mutations.fuzzy0.csv")
graph.prune(distance=3, sample_size=None, internal=True)
print(graph.number_sample_classes())
graph.to_tsv("/Users/rmnorris/Work/git/grapevine-global-tree/test/cog_global_2021-02-23_mutations.fuzzy3.tsv")
graph.output_sample_classes("/Users/rmnorris/Work/git/grapevine-global-tree/test/cog_global_2021-02-23_mutations.fuzzy3.csv")

