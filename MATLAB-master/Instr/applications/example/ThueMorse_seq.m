
function seq = ThueMorse_seq(seqlen)
    assert(seqlen >= 1, "sequence length must be larger than or equal to 1");
    seq = [0];
    while(length(seq) < seqlen)
        comp_seq = complement(seq);
        seq = cat(2, seq, comp_seq);
    end
    seq = seq(1:seqlen);
    function comp_l = complement(int_l)
        comp_l = [];
        for i = (1: length(int_l))
            if int_l(i) == 0
                comp_l(end+1) = 1;
            elseif int_l(i) == 1
                comp_l(end+1) = 0;
            end
        end
    end
end