using XSim
using JWAS
using DelimitedFiles
using CSV
using Statistics
using LinearAlgebra
using Plots
using Distributions
using Random
using Printf
pyplot()

function halfsibFamily(n)
    ind = 1:(1+2n)
    sire = [0;repeat([0],n);repeat([1],n)]
    dam  = [0;repeat([0],n);2:n+1]
    ped = [ind sire dam]
end

function extendPed(ped,ped1)
    n = ped[end,1]
    ind = ped1[:,1] .+ n
    sire = [s==0 ? 0 : s+n for s in ped1[:,2]]
    dam = [d==0 ? 0 : d+n for d in ped1[:,3]]
    ped = [ped;
           [ind sire dam]
        ]
end

function simData(nn,qtl,mrk,h2,founders)

    ped = halfsibFamily(nn[1])
    for i in nn[2:end]
        ped1 = halfsibFamily(i)
        ped = extendPed(ped,ped1)
    end

    writedlm("ped.csv",ped)

    pedigree   = get_pedigree("ped.csv",header=false,separator='\t');

    temp = Matrix(JWAS.PedModule.AInverse(pedigree));
    indx = [pedigree.idMap[string(i)].seqID for i in 1:size(ped,1)]
    Ai   = temp[indx,indx]
    A    = inv(Ai);

    animals = samplePed(ped,founders)
    QM = getOurGenotypes(animals)

    M = QM[:,mrk]
    Q = QM[:,qtl]
    n,p = size(M)
    return Q,M,A,ped
end

function getMME(nn,Q,M,A,ped,h2,wtMrk)
    
    n,p = size(M)
    
    sires = unique(ped[:,2])
    sires = sires[2:end]
    dams = unique(ped[:,3])
    dams = dams[2:end]
    offsp  = setdiff(ped[:,1],[sires;dams])
    cand = [1; sires]                       # candidates
    offsp = setdiff(offsp,cand)
    nOffsp = size(offsp,1)

    Z   = Matrix{Float64}(I,n,n)[offsp,:]
    ZPZ = Diagonal(Z'Z);

    G   = M*M'
    Het = mean(diag(G))
    G   = G/Het + I*0.001
    
    genVar = mean([dot(Q[i,:],Q[i,:]) for i=1:size(Q,1)])
    resVar = genVar*(1-h2)/h2

    var_g = genVar*wtMrk
    var_u = genVar*(1.0-wtMrk)

    H  = (G*var_g/genVar + A*var_u/genVar)
    genVarVec = diag(H)*genVar
    Hi = inv(H)
    J  = ones(nOffsp)

    λ  = resVar/genVar

    mme = [
    J'J J'Z
    Z'J ZPZ+Hi*λ   
    ]
    Ch = cholesky(Symmetric(mme));
    Pg = var_g/genVar*G*Hi
    Pu = var_u/genVar*A*Hi
    mmei = inv(Ch)
    d = diag(mmei)
    d = d[2:end]
    pev = d[cand]*resVar
    coruHatu = (genVarVec[cand] - pev) ./genVarVec[cand]
    coruHatu = map(x->x>0.0 ? sqrt(x) : 0.0,coruHatu)
    println("coruHatu")
    for i=1:size(nn,1)
        @printf("%4d %5.2f\n",nn[i],coruHatu[i])
    end
    return Pg,Pu,Ch,J,Z,offsp,cand,genVar,resVar,coruHatu
end

function simBLUP(Q,Pg,Pu,Ch,J,Z,n,offsp,cand,resVar)
    α = randn(size(Q,2))
    a = Q*α
    y = a + randn(n)*sqrt(resVar)
    y = y[offsp]
    rhs = [J'y; Z'y]
    sol = Ch\rhs
    aHat = sol[2:end]    # .+ sol[1]
    gHat = Pg*aHat
    uHat = Pu*aHat
    #a    = a    .- mean(a)
    #aHat = aHat .- mean(aHat)
    #gHat = gHat .- mean(gHat)
    #uHat = uHat .- mean(uHat)
    [aHat gHat uHat a][cand,:]
end

function accuracyBLUP(res,j)
    pred = vcat([i[j,:]' for i in res]...)
    cor_u_uHat = cor(pred)[1:3,4]
    covMat = cov(pred)
    d = diag(covMat)
    reg_u_uHat = covMat[1:3,4] ./ d[1:3]
    [cor_u_uHat; reg_u_uHat]
end

function varBLUP(res,j)
    pred = vcat([i[j,:]' for i in res]...)
    cor_u_uHat = cor(pred)[1:3,4]
    vars = var(pred,dims=1)
end

function covBLUP(res,j)
    pred = vcat([i[j,:]' for i in res]...)
    covs = cov(pred,dims=1)
    [covs[1:3,4]'./covs[4,4] covs[4,4]]
end

chrLength = 30.0
numChr    = 1
numLoci   = 6000
mutRate   = 0.0

mapPos = fill(0,60_000)
mrkPos = 10:20:60_000

qtlPos = mrkPos .+ 1

mapPos[mrkPos] .= 1
mapPos[qtlPos] .= 2

r = 30.00/60_000
map_pos = []
for (i,x) in enumerate(mapPos)
    if x>0
        push!(map_pos,i*r)
    end
end

mrkIndx = []
count = 0
for x in mapPos
    if x>0
        count += 1
        if x==1 
            push!(mrkIndx,count)
        end
    end
end

qtlIndx = []
count = 0
for x in mapPos
    if x>0
        count += 1
        if x==2 
            push!(qtlIndx,count)
        end
    end
end

geneFreq = fill(0.5,numLoci)
build_genome(numChr,chrLength,numLoci,float.(geneFreq),float.(map_pos),mutRate);

popSizeFounder = 100
sires = sampleFounders(popSizeFounder)
dams  = sampleFounders(popSizeFounder);

popSize = 200
nGen    = 100
sires,dams,gen = sampleRan(popSize, nGen,sires,dams;gen=1);

popSize = 1000
nGen    = 1
sires,dams,gen = sampleRan(popSize, nGen,sires,dams;gen=gen);

parents=concatCohorts(sires,dams);

Random.seed!(31415926123);
nn = [0,1,2,10,25,50,100,200,400]
h2 = 0.25
wtMrk = 1.0
Q,M,A,ped = simData(nn,qtlIndx,mrkIndx,h2,parents);

mean(Q)

Qc = Q .- mean(Q,dims=1)
Mc = M .- mean(M,dims=1);

wtMrk = 1.0
Pg,Pu,Ch,J,Z,offsp,cand,genVar,resVar,coruHatu =  getMME(nn,Qc,Mc,A,ped,h2,wtMrk);

nReps = 10_000
n   = size(Q,1)
res = [simBLUP(Qc,Pg,Pu,Ch,J,Z,n,offsp,cand,resVar) for i=1:nReps]
resAccuracy = [accuracyBLUP(res,j) for j=1:size(cand,1)]
corrReg = hcat(resAccuracy...)'
[nn corrReg]

plot(nn,[coruHatu corrReg[:,1]])

covMat = [nn vcat([covBLUP(res,i) for i=1:9]...)]

res100 = [coruHatu, corrReg, covMat];

wtMrk = 0.95
Pg,Pu,Ch,J,Z,offsp,cand,genVar,resVar,coruHatu =  getMME(nn,Qc,Mc,A,ped,h2,wtMrk);

nReps = 10_000
n   = size(Q,1)
res = [simBLUP(Qc,Pg,Pu,Ch,J,Z,n,offsp,cand,resVar) for i=1:nReps]
resAccuracy = [accuracyBLUP(res,j) for j=1:size(cand,1)]
corrReg = hcat(resAccuracy...)'
[nn corrReg]

plot(nn,[coruHatu corrReg[:,1]])

covMat = [nn vcat([covBLUP(res,i) for i=1:9]...)]

res95 = [coruHatu, corrReg, covMat];

wtMrk = 0.5
Pg,Pu,Ch,J,Z,offsp,cand,genVar,resVar,coruHatu =  getMME(nn,Qc,Mc,A,ped,h2,wtMrk);

nReps = 10_000
n   = size(Q,1)
res = [simBLUP(Qc,Pg,Pu,Ch,J,Z,n,offsp,cand,resVar) for i=1:nReps]
resAccuracy = [accuracyBLUP(res,j) for j=1:size(cand,1)]
corrReg = hcat(resAccuracy...)'
[nn corrReg]

plot(nn,[coruHatu corrReg[:,1]])

covMat = [nn vcat([covBLUP(res,i) for i=1:9]...)]

res50 = [coruHatu, corrReg, covMat];

wtMrk = 0.0
Pg,Pu,Ch,J,Z,offsp,cand,genVar,resVar,coruHatu =  getMME(nn,Qc,Mc,A,ped,h2,wtMrk);

nReps = 10_000
n   = size(Q,1)
res = [simBLUP(Qc,Pg,Pu,Ch,J,Z,n,offsp,cand,resVar) for i=1:nReps]
resAccuracy = [accuracyBLUP(res,j) for j=1:size(cand,1)]
corrReg = hcat(resAccuracy...)'
[nn corrReg]

plot(nn,[coruHatu corrReg[:,1]])

covMat = [nn vcat([covBLUP(res,i) for i=1:9]...)]
