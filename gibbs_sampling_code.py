import random

#Random kmer
def random_kmer(seq,k):
    start=random.randint(0,len(seq)-k)
        
    return seq[start:start+k]

#Build_profile
def build_profile(motif):
    profile=[]
    n=len(motif[0])
    for i in range(n):
      count={'A':1,'T':1,'G':1,'C':1}
      for m in motif:
          count[m[i]]+=1
      total=sum(count.values())
      prob={}
      for x in 'ATGC':
        prob[x]=count[x]/total
      profile.append(prob)
    return profile

#calculate_probability
def get_prob(kmer,profile):
    p=1
    for i in range(len(kmer)):
        ch=kmer[i]
        p=p*profile[i][ch]
    return p

#Profile_random_kmer
def profile_random_kmer(seq,k,profile):
    kmers=[]
    probs=[]
    for i in range(len(seq)-k+1):
        kmer=seq[i:i+k]
        kmers.append(kmer)
        p=get_prob(kmer,profile)
        probs.append(p)
    selected=random.choices(kmers,weights=probs,k=1)

    return selected[0]

# Score function
def find_score(motif):
    score=0
    n=len(motif)
    k=len(motif[0])
    for i in range(k):
        A=T=G=C=0
        for m in motif:
            if(m[i]=='A'):
                A+=1
            if(m[i]=='T'):
                T+=1
            if(m[i]=='G'):
                G+=1
            if(m[i]=='C'):
                C+=1
        max_count=max(A,T,G,C)
        score+=(n-max_count)
    
    return score

#Gibbs function
def Gibbs(DNA,K,N=1000):
    motif=[]
    for seq in DNA:
        motif.append(random_kmer(seq,K))
    
    best_motif=motif[:]

    for step in range(N):

        i=random.randint(0,len(DNA)-1)

        temp=motif[:]
        temp.pop(i)

        profile=build_profile(temp)

        new_motif=profile_random_kmer(DNA[i],k,profile)
        motif[i]=new_motif

        if find_score(motif)<find_score(best_motif):
            best_motif=motif[:]

    return best_motif

n=int(input("Enter the no of DNA sequence:"))
DNA=[]
for i in range(n):
    seq=input("Enter DNA Squence:")
    DNA.append(seq)
k=int(input("Enter K value:"))
ans=Gibbs(DNA,k)
print("The motif set is:",ans)
print("The Score is:",find_score(ans))