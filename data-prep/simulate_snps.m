function output_table=simulate_snps(dependency_directory,systematic_name,dna_sequence)

%some sequences seem to not be mod(length,3)==0
n_bases=3*floor(length(dna_sequence)/3)-3; %don't do stop codon

base_array={'A','C','G','T'};

initial_base=cell(3*n_bases,1);
new_base=initial_base;

initial_residue=initial_base;
new_residue=initial_base;

initial_codon=initial_base;
new_codon=initial_base;

is_mis=zeros(3*n_bases,1);
residue_number=is_mis;

is_ts=is_mis;

m=1;
for i=1:n_bases

    temp_base=dna_sequence(i);
    residue_idx=(1:3)+(ceil(i/3)-1)*3;
    temp_residue=nt2aa(dna_sequence(residue_idx),'AlternativeStartCodons',false);
    
    mutation_array=base_array(~ismember(base_array,temp_base));
    for j=1:length(mutation_array)
        
        initial_base{m}=temp_base;
        new_base{m}=mutation_array{j};
        
        new_sequence=dna_sequence;
        new_sequence(i)=mutation_array{j};
        
        initial_residue{m}=temp_residue;
        new_residue{m}=nt2aa(new_sequence(residue_idx),'AlternativeStartCodons',false);
        
        initial_codon{m}=dna_sequence(residue_idx);
        new_codon{m}=new_sequence(residue_idx);
        
        is_mis(m)=~strcmp(initial_residue{m},new_residue{m});
        
        %Ts/Tv
        if (strcmp(initial_base{m},'A')&&strcmp(new_base{m},'G'))||...
                (strcmp(initial_base{m},'G')&&strcmp(new_base{m},'A'))||...
                (strcmp(initial_base{m},'C')&&strcmp(new_base{m},'T'))||...
                (strcmp(initial_base{m},'T')&&strcmp(new_base{m},'C'))
            is_ts(m)=1;
        end
        
        residue_number(m)=ceil(i/3);
        
        m=m+1;
        
    end
    
end

output_table=table(residue_number,initial_base,new_base,initial_residue,new_residue,...
    initial_codon,new_codon,is_mis,is_ts);

save([dependency_directory 'mutation-tables/' systematic_name '_mutation_table.mat'],'output_table')

end

