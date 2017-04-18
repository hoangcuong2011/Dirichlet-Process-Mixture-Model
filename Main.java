/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
import cc.mallet.types.Dirichlet;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
/**
 *
 * @author hoangcuong2011
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    
    static int burnin = 10;
    
    static Random r = new Random();
    static int[] Class;
    static double[] Samples;
    
    static int number_of_clusters = 5;
    static double varianceFIX = 0.001;
    
    private static double getGaussian(double aMean, double aVariance) {
        
        return aMean + r.nextGaussian() * aVariance;
    }

    public static double[] generateSamplesforToyProblems(int number_of_clusters) {
        int points_per_cluster = 50;
        Class = new int[points_per_cluster*number_of_clusters];
        double[] samples = new double[points_per_cluster*number_of_clusters];
        int count = 0;
        for (int cluster = 0; cluster < number_of_clusters; cluster++) {
            for (int i = 0; i < points_per_cluster; i++) {
                samples[count] = getGaussian((double) (cluster+1)/(double) number_of_clusters, Math.sqrt(varianceFIX));
                Class[count] = cluster;
                count++;
            }
        }
        return samples;
    }

    public static double probabilitydensityfunct(double mean, double variance, double y) {
        double result = 1.0 / Math.sqrt(2.0 * Math.PI * variance) * Math.exp(-1 * (y - mean) * (y - mean) / 2 / variance);
        return result;
    }

    public static int[] findRank(double[] inp) {
        int[] outp = new int[inp.length];
        for (int i = 0; i < inp.length; i++) {
            for (int k = 0; k < inp.length; k++) {
                if (inp[k] < inp[i]) {
                    outp[i]++;
                }
            }
        }
        return outp;
    }

    public static void EM(int number_of_clusters, double[] Samples, int[] Class) {
        double likelihoodMAX = -1000000;
        double accuracy = -1;

        for (int trial = 0; trial < 1000; trial++) {
            double lambdas[] = new double[number_of_clusters];
            double means[] = new double[number_of_clusters];
            double variance[] = new double[number_of_clusters];
            for (int i = 0; i < number_of_clusters; i++) {
                lambdas[i] = 1.0 / number_of_clusters;
                means[i] = r.nextDouble();

                variance[i] = varianceFIX;

            }

            for (int iteration = 0; iteration < 8; iteration++) {
                double cache_mean[] = new double[number_of_clusters];
                double total_statistics[] = new double[number_of_clusters];
                double cache_variance[] = new double[number_of_clusters];
                double cache_lambdas[] = new double[number_of_clusters];
                //Estep
                double likelihood = 0.0;
                double truecase = 0;
                int[] Ranks = findRank(means);
                int myclass[] = new int[Samples.length];
                for (int j = 0; j < Samples.length; j++) {
                    double r_j[] = new double[number_of_clusters];
                    double sum = 0.0;
                    double dmax = -1.0;
                    int imax = -1;
                    for (int i = 0; i < number_of_clusters; i++) {
                        r_j[i] = lambdas[i] * probabilitydensityfunct(means[i], variance[i], Samples[j]);
                        if (r_j[i] > dmax) {
                            dmax = r_j[i];
                            imax = i;
                        }
                        if (r_j[i] < 1e-200) {
                            r_j[i] = 1e-200;
                        }
                        if (r_j[i] <= 0) {
                            System.out.println("Shit happens");
                        }
                        sum += r_j[i];
                    }
                    myclass[j] = imax;

                    if (Ranks[imax] == Class[j]) {
                        truecase++;
                    }
                    if (sum < 1e-200) {
                        sum = 1e-200;
                    }
                    likelihood += Math.log(sum);
                    for (int i = 0; i < number_of_clusters; i++) {
                        r_j[i] /= sum;
                    }

                    //collecting statistics
                    for (int i = 0; i < number_of_clusters; i++) {
                        cache_mean[i] += r_j[i] * Samples[j];
                        total_statistics[i] += r_j[i];
                        cache_variance[i] += r_j[i] * (Samples[j] - means[i]) * (Samples[j] - means[i]);
                        cache_lambdas[i] += r_j[i];
                    }
                }
                //Mstep
                for (int i = 0; i < number_of_clusters; i++) {
                    means[i] = cache_mean[i] / total_statistics[i];
                    //variance[i] = cache_variance[i] / total_statistics[i];
                    lambdas[i] = cache_lambdas[i] / Samples.length;
                }
                if (iteration % 1 == 0) {
                    //System.out.println("\nIteration " + iteration);
                    //System.out.println("Log-Likelihood ..." + likelihood);
                    truecase = measure(Class, myclass);
                    //System.out.println("Accuracy ..." + ((double) (truecase)));
                    //System.out.println("means ...");
                    //for (int i = 0; i < number_of_clusters; i++) {
                        //System.out.print(" " + means[i]);
                    //}
                    //System.out.println();
                    if (likelihood > likelihoodMAX) {
                        likelihoodMAX = likelihood;
                        accuracy = ((double) (truecase));
                    }
                    /*System.out.println("variances ...");
                for (int i = 0; i < number_of_clusters; i++) {
                    System.out.print(" " + variance[i]);
                }
                System.out.println();*/
                    //System.out.println("lambdas ...");
                    //for (int i = 0; i < number_of_clusters; i++) {
                        //System.out.print(" " + lambdas[i]);
                    //}
                }
            }

        }
        System.out.println("\nACCURACY " + accuracy);
    }

    public static int[] GibbsdiscreteSamples(int number_of_cluster, double[] Samples, double[] thetas, 
            double[] variances, int[] classCountings, int[] originalSequencesZ) {
        int[] sequencesZ = new int[Samples.length];
        double As[][] = new double[Samples.length][];
        Random r = new Random();
        for (int sample = 0; sample < Samples.length; sample++) {
            As[sample] = new double[number_of_cluster];
            int k = originalSequencesZ[sample];
            classCountings[k]-=1;
            double[] Probs = new double[number_of_cluster];
            for (int i = 0; i < number_of_cluster; i++) {
                Probs[i] = probabilitydensityfunct(thetas[i], variances[i], Samples[sample]);
                Probs[i] *= (classCountings[i]+1.0/(double) number_of_cluster)/(Samples.length+1.0-1.0);
            }
            //normalizeProbs first
            double sum = 0;
            for (int i = 0; i < number_of_cluster; i++) {
                sum += Probs[i];
            }
            for (int i = 0; i < number_of_cluster; i++) {
                Probs[i] /= sum;
            }
            for (int i = 0; i < number_of_cluster; i++) {
                if (i > 0) {
                    As[sample][i] = As[sample][i - 1] + Probs[i];
                } else {
                    As[sample][i] = Probs[i];
                }
            }
            double rndNumber = r.nextDouble();
            LOOP:
            for (int i = 0; i < number_of_cluster; i++) {
                if (i == 0) {
                    if (rndNumber <= As[sample][i]) {
                        //pick i
                        sequencesZ[sample] = i;
                        classCountings[i]+=1;
                        break LOOP;
                    }
                } else {
                    if (rndNumber > As[sample][i - 1] && rndNumber <= As[sample][i]) {
                        sequencesZ[sample] = i;
                        classCountings[i]+=1;
                        break LOOP;
                    }
                }
            }
        }
        
        return sequencesZ;
    }
    
    
    
    
    

    public static  int[] discreteSamples(int number_of_cluster, double[] Samples, double[] thetas, double[] pis, double[] variances) {
        
        int[] sequencesZ = new int[Samples.length];
        double As[][] = new double[Samples.length][];
        for (int sample = 0; sample < Samples.length; sample++) {
            As[sample] = new double[number_of_cluster];
            double[] Probs = new double[number_of_cluster];
            for (int i = 0; i < number_of_cluster; i++) {
                Probs[i] = probabilitydensityfunct(thetas[i], variances[i], Samples[sample]);
                Probs[i] *= pis[i];
            }
            //normalizeProbs first
            double sum = 0;
            for (int i = 0; i < number_of_cluster; i++) {
                sum+= Probs[i];
            }
            for (int i = 0; i < number_of_cluster; i++) {
                Probs[i]/=sum;
            }
            for (int i = 0; i < number_of_cluster; i++) {
                if (i > 0) {
                    As[sample][i] = As[sample][i - 1] + Probs[i];
                } else {
                    As[sample][i] = Probs[i];
                }
            }
        }
        Random r = new Random();
        for (int sample = 0; sample < Samples.length; sample++) {            
            double rndNumber = r.nextDouble();
            LOOP: for (int i = 0; i < number_of_cluster; i++) {
                if (i == 0) {
                    if (rndNumber <= As[sample][i]) {
                        //pick i
                        sequencesZ[sample] = i;
                        break LOOP;
                    }
                }
                else {
                    if (rndNumber > As[sample][i-1]&&rndNumber<=As[sample][i]) {
                        sequencesZ[sample] = i;
                        break LOOP;
                    }
                }
            }
        }
        return sequencesZ;
    }
    public static void GibbsSamplingforFMM(int number_of_clusters, double[] Samples, int[] Class) {
        System.out.println("GibbsSamplingforFMM");
        double pis[] = new double[number_of_clusters];
        double thetas[] = new double[number_of_clusters];
        double variances[] = new double[number_of_clusters];
        double cacheZs[][] = new double[Samples.length][];
        double cacheThetas[] = new double[number_of_clusters];
        double cachePis[] = new double[number_of_clusters];
        for(int i = 0; i < Samples.length; i++)
        {
            cacheZs[i] = new double[number_of_clusters];
        }
        for(int i = 0; i < number_of_clusters; i++)
        {
            variances[i] = varianceFIX;
        }
        Random r = new Random();
        ArrayList<ArrayList<Double>> list = new ArrayList<ArrayList<Double>>();
        for(int i = 0; i < number_of_clusters; i++) {
            pis[i] = 1.0/number_of_clusters;        
            list.add(new ArrayList<>());
        }
        
        for(int i = 0; i < number_of_clusters; i++)
            thetas[i] = r.nextDouble();
        int counts[] = new int[number_of_clusters];
        
        double maximumiteration = 3000;
        for (int iteration = 0; iteration < maximumiteration; iteration++) {
            //sampling Z
            int[] sequencesZ = discreteSamples(number_of_clusters, Samples, thetas, pis, variances);
            if (iteration > burnin) {
                if ((iteration - burnin) % 100==0) {                    
                    for (int i = 0; i < sequencesZ.length; i++) {
                        cacheZs[i][sequencesZ[i]] = cacheZs[i][sequencesZ[i]] + 1;
                    }
                    for (int i = 0; i < cacheThetas.length; i++) {
                        counts[i]+=1;
                        cacheThetas[i] += thetas[i];
                        list.get(i).add(thetas[i]);
                        cachePis[i] +=pis[i];
                    }
                }
            }
            if(iteration==(maximumiteration-1)) {
                
                double truecase = 0;
                int[] myclass = new int[Samples.length];
                for (int i = 0; i < sequencesZ.length; i++) {
                    
                    double dmax = -1.0;
                    int imax = -1;
                    for (int j = 0; j < number_of_clusters; j++) {
                        
                        /*double tempvalue = probabilitydensityfunct(cacheThetas[j]/(double) counts[j], varianceFIX, Samples[i]);
                        tempvalue *= (cachePis[j]/(double) counts[j]);*/
                        
                        double tempvalue = cacheZs[i][j];
                        if (dmax < tempvalue) {
                            dmax = tempvalue;
                            imax = j;
                        }                        
                    }
                    myclass[i] = imax;
                }
                
                
                truecase = measure(Class, myclass);
                System.out.println("Iteration "+iteration+"\nAccuracy ..."+((double) (truecase)));
                /*System.out.println("means ...");
                for (int i = 0; i < number_of_clusters; i++) {
                    System.out.print(" " + cacheThetas[i]/(double) counts[i]+"~"+cachePis[i]/(double) counts[i]);
                }*/
                System.out.println();
                if (list.get(0).size() > 100) {
                    double subsets[] = new double[100];
                    for (int i = 0; i < 100; i++) {
                        subsets[i] = list.get(0).get(list.get(0).size() - i - 1);//pick cluster 0 as sanity check
                        //System.out.println(subsets[i]);
                    }
                    Statistics mypackage = new Statistics(subsets);
                    //System.out.println("variance " + mypackage.getVariance());
                } else {
                    double subsets[] = new double[list.get(0).size()];
                    for (int i = 0; i < list.get(0).size(); i++) {
                        subsets[i] = list.get(0).get(list.get(0).size() - i - 1);//pick cluster 0 as sanity check
                        //System.out.println(subsets[i]);
                    }
                    Statistics mypackage = new Statistics(subsets);
                    //System.out.println("variance " + mypackage.getVariance());
                }

                /*System.out.println("means ...");
                for (int i = 0; i < number_of_clusters; i++) {
                    System.out.print(" " + thetas[i]);
                }
                System.out.println();
                System.out.println("variances ...");
                for (int i = 0; i < number_of_clusters; i++) {
                    System.out.print(" " + variances[i]);
                }
                System.out.println();
                System.out.println("lambdas ...");
                for (int i = 0; i < number_of_clusters; i++) {
                    System.out.print(" " + pis[i]);
                }*/
            }
            double[] cacheStatistics = new double[number_of_clusters];
            ArrayList<ArrayList<Double>> clusteringPoints = new ArrayList<ArrayList<Double>>();
            for (int i = 0; i < number_of_clusters; i++) {
                cacheStatistics[i] = 1.0/(double) number_of_clusters;//alpha = 1
                clusteringPoints.add(new ArrayList<>());
            }
            for (int i = 0; i < sequencesZ.length; i++) {
                cacheStatistics[sequencesZ[i]] = cacheStatistics[sequencesZ[i]] + 1;
                clusteringPoints.get(sequencesZ[i]).add(Samples[i]);
            }
            for (int i = 0; i < number_of_clusters; i++) {
                if(cacheStatistics[i]==1.0/(double)number_of_clusters) {
                    //System.out.println("Shit happens ....");
                }
            }
            Dirichlet Dirichletsampling = new Dirichlet(cacheStatistics);
            pis = Dirichletsampling.nextDistribution();

            double newMeansforsampling[] = new double[number_of_clusters];
            double newVarianceforsampling[] = new double[number_of_clusters];
            double mu_0 = 0.0;
            double variance_0 = 1000;
            for (int i = 0; i < number_of_clusters; i++) {
                /*another formula
                newMeansforsampling[i] = variance_0 * means[i] / (varianceFIX/ ((double) clusteringPoints.get(i).size()+0.0001) + variance_0);
                newMeansforsampling[i] = newMeansforsampling[i]
                        + varianceFIX/ (varianceFIX / (double) clusteringPoints.get(i).size() +0.0001 + variance_0) * mu_0; 
                */
                double sum = 0;
                for (int j = 0; j < clusteringPoints.get(i).size(); j++) {
                    sum += clusteringPoints.get(i).get(j);
                }
                
                newMeansforsampling[i] = varianceFIX*mu_0+variance_0*sum;
                newMeansforsampling[i] /= (varianceFIX+variance_0*clusteringPoints.get(i).size()); 
                
                //I ignore this because I assume mu_0 = 0 
                newVarianceforsampling[i] = varianceFIX*variance_0;
                newVarianceforsampling[i] /= (varianceFIX+variance_0*clusteringPoints.get(i).size());
            }
            
            for (int i = 0; i < number_of_clusters; i++) {
                if(Double.isInfinite(newMeansforsampling[i]) || Double.isNaN(newMeansforsampling[i])) {
                    System.out.println("Shit happens ....");
                    //thetas[i] =0.001+r.nextDouble()/100.0;
                    continue;
                }
                thetas[i] = getGaussian(newMeansforsampling[i], Math.sqrt(newVarianceforsampling[i]));
            }
        }

    }
    
    
    
    public static void CollapsedGibbsSamplingforFMM(int number_of_clusters, double[] Samples, int[] Class) {
        System.out.println("CollapsedGibbsSamplingforFMM");
        //double pis[] = new double[number_of_clusters];
        double thetas[] = new double[number_of_clusters];
        double variances[] = new double[number_of_clusters];
        double cacheZs[][] = new double[Samples.length][];
        double cacheThetas[] = new double[number_of_clusters];
        double cachePis[] = new double[number_of_clusters];
        for(int i = 0; i < Samples.length; i++)
        {
            cacheZs[i] = new double[number_of_clusters];
        }
        for(int i = 0; i < number_of_clusters; i++)
        {
            variances[i] = varianceFIX;
        }
        Random r = new Random();
        ArrayList<ArrayList<Double>> list = new ArrayList<ArrayList<Double>>();
        for(int i = 0; i < number_of_clusters; i++) {
            //pis[i] = 1.0/number_of_clusters;        
            list.add(new ArrayList<>());
        }
        
        for(int i = 0; i < number_of_clusters; i++)
            thetas[i] = r.nextDouble();
        int counts[] = new int[number_of_clusters];
        
        int classCountings[] = new int[number_of_clusters];
        int[] sequencesZ = new int[Samples.length];
        Random mr = new Random();
        for(int i = 0; i < sequencesZ.length; i++) {
            int classVariable = mr.nextInt(number_of_clusters);
            sequencesZ[i] = classVariable;
            classCountings[classVariable]+=1;
        }
        double maximumiteration = 3000;
        for (int iteration = 0; iteration < maximumiteration; iteration++) {
            //sampling Z
            sequencesZ = GibbsdiscreteSamples(number_of_clusters, Samples, thetas, variances, classCountings, sequencesZ);
            if (iteration > burnin) {
                if ((iteration - burnin) % 100==0) {                    
                    for (int i = 0; i < sequencesZ.length; i++) {
                        cacheZs[i][sequencesZ[i]] = cacheZs[i][sequencesZ[i]] + 1;
                    }
                    for (int i = 0; i < cacheThetas.length; i++) {
                        counts[i]+=1;
                        cacheThetas[i] += thetas[i];
                        list.get(i).add(thetas[i]);
                        //cachePis[i] +=pis[i];
                    }
                }
            }
            if(iteration==(maximumiteration-1)) {
                
                double truecase = 0;
                int[] myclass = new int[Samples.length];
                for (int i = 0; i < sequencesZ.length; i++) {
                    
                    double dmax = -1.0;
                    int imax = -1;
                    for (int j = 0; j < number_of_clusters; j++) {
                        
                        /*double tempvalue = probabilitydensityfunct(cacheThetas[j]/(double) counts[j], varianceFIX, Samples[i]);
                        tempvalue *= (cachePis[j]/(double) counts[j]);*/
                        
                        double tempvalue = cacheZs[i][j];
                        if (dmax < tempvalue) {
                            dmax = tempvalue;
                            imax = j;
                        }                        
                    }
                    myclass[i] = imax;
                }
                
                
                truecase = measure(Class, myclass);
                System.out.println("Iteration "+iteration+"\nAccuracy ..."+((double) (truecase)));
                System.out.println("means ...");
                for (int i = 0; i < number_of_clusters; i++) {
                    System.out.print(" " + cacheThetas[i]/(double) counts[i]+"~"+cachePis[i]/(double) counts[i]);
                }
                System.out.println();
                if (list.get(0).size() > 100) {
                    double subsets[] = new double[100];
                    for (int i = 0; i < 100; i++) {
                        subsets[i] = list.get(0).get(list.get(0).size() - i - 1);//pick cluster 0 as sanity check
                        //System.out.println(subsets[i]);
                    }
                    Statistics mypackage = new Statistics(subsets);
                    //System.out.println("variance " + mypackage.getVariance());
                } else {
                    double subsets[] = new double[list.get(0).size()];
                    for (int i = 0; i < list.get(0).size(); i++) {
                        subsets[i] = list.get(0).get(list.get(0).size() - i - 1);//pick cluster 0 as sanity check
                        //System.out.println(subsets[i]);
                    }
                    Statistics mypackage = new Statistics(subsets);
                    //System.out.println("variance " + mypackage.getVariance());
                }

                /*System.out.println("means ...");
                for (int i = 0; i < number_of_clusters; i++) {
                    System.out.print(" " + thetas[i]);
                }
                System.out.println();
                System.out.println("variances ...");
                for (int i = 0; i < number_of_clusters; i++) {
                    System.out.print(" " + variances[i]);
                }
                System.out.println();
                System.out.println("lambdas ...");
                for (int i = 0; i < number_of_clusters; i++) {
                    System.out.print(" " + pis[i]);
                }*/
            }
            double[] cacheStatistics = new double[number_of_clusters];
            ArrayList<ArrayList<Double>> clusteringPoints = new ArrayList<ArrayList<Double>>();
            for (int i = 0; i < number_of_clusters; i++) {
                cacheStatistics[i] = 1.0/(double) number_of_clusters;//alpha = 1
                clusteringPoints.add(new ArrayList<>());
            }
            for (int i = 0; i < sequencesZ.length; i++) {
                cacheStatistics[sequencesZ[i]] = cacheStatistics[sequencesZ[i]] + 1;
                clusteringPoints.get(sequencesZ[i]).add(Samples[i]);
            }
            for (int i = 0; i < number_of_clusters; i++) {
                if(cacheStatistics[i]==1.0/(double)number_of_clusters) {
                    //System.out.println("Shit happens ....");
                }
            }
            //Dirichlet Dirichletsampling = new Dirichlet(cacheStatistics);
            //pis = Dirichletsampling.nextDistribution();

            double newMeansforsampling[] = new double[number_of_clusters];
            double newVarianceforsampling[] = new double[number_of_clusters];
            double mu_0 = 0.0;
            double variance_0 = 1000;
            for (int i = 0; i < number_of_clusters; i++) {
                /*another formula
                newMeansforsampling[i] = variance_0 * means[i] / (varianceFIX/ ((double) clusteringPoints.get(i).size()+0.0001) + variance_0);
                newMeansforsampling[i] = newMeansforsampling[i]
                        + varianceFIX/ (varianceFIX / (double) clusteringPoints.get(i).size() +0.0001 + variance_0) * mu_0; 
                */
                double sum = 0;
                for (int j = 0; j < clusteringPoints.get(i).size(); j++) {
                    sum += clusteringPoints.get(i).get(j);
                }
                
                newMeansforsampling[i] = varianceFIX*mu_0+variance_0*sum;
                newMeansforsampling[i] /= (varianceFIX+variance_0*clusteringPoints.get(i).size()); 
                
                //I ignore this because I assume mu_0 = 0 
                newVarianceforsampling[i] = varianceFIX*variance_0;
                newVarianceforsampling[i] /= (varianceFIX+variance_0*clusteringPoints.get(i).size());
            }
            
            for (int i = 0; i < number_of_clusters; i++) {
                if(Double.isInfinite(newMeansforsampling[i]) || Double.isNaN(newMeansforsampling[i])) {
                    System.out.println("Shit happens ....");
                    //thetas[i] =0.001+r.nextDouble()/100.0;
                    continue;
                }
                thetas[i] = getGaussian(newMeansforsampling[i], Math.sqrt(newVarianceforsampling[i]));
            }
        }

    }

    static class ClusterDataStructure {
        int index = 0;//index - name of cluster
        boolean[] listID; //I save the element id instead for simplicity
        double theta = 0; //mean of the cluster
        int count = 0; //number of elements in the cluster
        public ClusterDataStructure(int SampleSize) {
            listID = new boolean[SampleSize];//all of them are false            
        }
    }


    
    public static int[] GibbsdiscreteSamplesDPMM(double[] Samples, ArrayList<ClusterDataStructure> list, 
            int[] originalSequencesZ, int totalLimitClustersToDate) {
        int[] sequencesZ = new int[Samples.length];
        Random r = new Random();
        for (int sample = 0; sample < Samples.length; sample++) {
            double As[] = new double[totalLimitClustersToDate];
            int k = originalSequencesZ[sample];
            list.get(k).listID[sample] = false; //I remove it to do sample a new cluster.
            list.get(k).count = list.get(k).count - 1;
            if(list.get(k).count==0) {
                //ok I remove it
                list.set(k, null);
            }
            double[] Probs = new double[totalLimitClustersToDate];
            boolean nullforthefirstime = false;
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                if (list.get(i) != null) {
                    Probs[i] = probabilitydensityfunct(list.get(i).theta, varianceFIX, Samples[sample]);
                    Probs[i] *= (list.get(i).count) / (Samples.length + 1.0 - 1.0);
                } else {
                    //if null for the first time; ok new cluster
                    if(nullforthefirstime==false) {
                        nullforthefirstime = true; //i do this onlly once.
                        Probs[i] = integralFormula(Samples[sample], varianceFIX, 0, 100);
                        Probs[i] *= 1.0 / (Samples.length + 1.0 - 1.0);
                    }
                }                
            }
            //normalizeProbs first
            double sum = 0;
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                sum += Probs[i];
            }
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                Probs[i] /= sum;
            }
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                if (i > 0) {
                    As[i] = As[i - 1] + Probs[i];
                } else {
                    As[i] = Probs[i];
                }
            }
            double rndNumber = r.nextDouble();
            LOOP:
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                if (i == 0) {
                    if (rndNumber <= As[i]) {
                        //pick i
                        sequencesZ[sample] = i;
                        if(list.get(i)==null) {
                            list.set(i, new ClusterDataStructure(Samples.length));
                            list.get(i).index = i;
                            //sample new mean - refer to jordan - need to proof

                            double newmean = varianceFIX * 0 + 100 * Samples[sample];
                            newmean /= (varianceFIX + 100);

                            //I ignore this because I assume mu_0 = 0 
                            double newvariance = varianceFIX * 100;
                            newvariance /= (varianceFIX + 100);
                            list.get(i).theta = getGaussian(newmean, Math.sqrt(newvariance));
                        }
                        list.get(i).listID[sample] = true;
                        list.get(i).count = list.get(i).count +1;
                        break LOOP;
                    }
                } else {
                    if (rndNumber > As[i - 1] && rndNumber <= As[i]) {            
                        if(i>25) {
                            System.out.println("Probably shit happens ...");
                        }
                        sequencesZ[sample] = i;
                        if(list.get(i)==null) {
                            list.set(i, new ClusterDataStructure(Samples.length));
                            list.get(i).index = i;
                            
                            double newmean = varianceFIX * 0 + 100 * Samples[sample];
                            newmean /= (varianceFIX + 100);

                            //sample new mean - refer to jordan - need to proof
                            double newvariance = varianceFIX * 100;
                            newvariance /= (varianceFIX + 100);
                            list.get(i).theta = getGaussian(newmean, Math.sqrt(newvariance));
                        }
                        list.get(i).listID[sample] = true;
                        list.get(i).count = list.get(i).count +1;
                        break LOOP;
                    }
                }
            }
        }
        return sequencesZ;
    }
    
    public static void GibbsSamplingforDGMM(double[] Samples, int[] Class) {
        System.out.println("GibbsSamplingforDGMM");
        double sumStatistics = 0;
        int countTimes = 0;
        ArrayList<ArrayList<Integer>> cacheZs = new ArrayList<>(); //an arraylist of arraylist with sample id and the list of all assined cluster
        for(int i = 0; i < Samples.length; i++)
        {
            cacheZs.add(new ArrayList<>());
        }
        
        Random r = new Random();
        ArrayList<ClusterDataStructure> list = new ArrayList<>();
        int maximumcluster = 1000;//just a stupid number
        for(int i = 0; i < maximumcluster; i++) {
            list.add(null);
        }
        list.set(0, new ClusterDataStructure(Samples.length));
        list.get(0).index = 0;        
        Random mr = new Random();
        double sum = 0;
        for(int i = 0; i < Samples.length; i++) {
            list.get(0).listID[i] = true;
            sum+=Samples[i];
        }
        list.get(0).count = Samples.length;
        //list.get(0).theta = r.nextDouble(); //in principle we don't need to specify this?!
        list.get(0).theta = sum/Samples.length; //in principle we don't need to specify this?!
        
        
        int[] sequencesZ = new int[Samples.length];
        //at the begining, all sequences Z sample to cluster 0
        
        
        double maximumiteration = 3000;
        for (int iteration = 0; iteration < maximumiteration; iteration++) {
            //sampling Z
            sequencesZ = GibbsdiscreteSamplesDPMM(Samples, list, sequencesZ, maximumcluster);
            if (iteration > burnin) {
                if ((iteration - burnin) % 1==0) {                    
                    for (int i = 0; i < sequencesZ.length; i++) {
                        cacheZs.get(i).add(sequencesZ[i]);
                    }
                }
            
            if(iteration==(maximumiteration-1) || (iteration - burnin) % 100==0) {         
                double truecase = 0;
                int[] myclass = new int[Samples.length];
                for (int i = 0; i < sequencesZ.length; i++) {
                    
                    int[] counts = new int[maximumcluster];
                    if (cacheZs.get(i).size() > 100) { //not a great solution, I know
                        for (int j = cacheZs.get(i).size()-100; j < cacheZs.get(i).size(); j++) {
                            counts[cacheZs.get(i).get(j)] = counts[cacheZs.get(i).get(j)] + 1;
                        }
                    } else {
                        for (int j = 0; j < cacheZs.get(i).size(); j++) {
                            counts[cacheZs.get(i).get(j)] = counts[cacheZs.get(i).get(j)] + 1;
                        }
                    }
                    double dmax = -1.0;
                    int imax = -1;
                    for (int j = 0; j < maximumcluster; j++) {
                        double tempvalue = counts[j];
                        if(i == 100 && (double)cacheZs.get(i).size()>10 && tempvalue>0) {
                            //System.out.println("abc "+j+"~"+tempvalue/(double)cacheZs.get(i).size()*100);
                        }
                        if (dmax < tempvalue) {
                            dmax = tempvalue;
                            imax = j;
                        }
                    }
                    myclass[i] = imax;
                }
                truecase = measure(Class, myclass);
                
                countTimes++;
                sumStatistics+=truecase;
                
                //System.out.println("Iteration " + iteration + "\nAccuracy ..." + ((double) (truecase)));
                int activeclass = 0;
                for (int i = 0; i < maximumcluster; i++) {
                    if (list.get(i) == null) {
                        continue;
                    }
                    activeclass++;

                }
                System.out.println("Iteration " + iteration + " Active classes" + activeclass);
                
                
                for (int i = 0; i < maximumcluster; i++) {
                    if (list.get(i) == null) {
                        continue;
                    }
                    System.out.print("~" + i);

                }
                System.out.println("");

                }
                double mu_0 = 0.0;
                double variance_0 = 1000;

                for (int i = 0; i < maximumcluster; i++) {
                    if (list.get(i) == null) {
                        continue;
                    }
                    sum = 0;
                    for (int j = 0; j < list.get(i).listID.length; j++) {
                        if (list.get(i).listID[j] == true) {
                            sum += Samples[j];
                        }
                    }

                    double newMeansforsampling = varianceFIX * mu_0 + variance_0 * sum;
                    newMeansforsampling /= (varianceFIX + variance_0 * list.get(i).count);

                    //I ignore this because I assume mu_0 = 0 
                    double newVarianceforsampling = varianceFIX * variance_0;
                    newVarianceforsampling /= (varianceFIX + variance_0 * list.get(i).count);
                    list.get(i).theta = getGaussian(newMeansforsampling, Math.sqrt(newVarianceforsampling));
                }
                //System.out.println("Iteration "+iteration+" Active classes" + activeclass);

            }

        }
        System.out.println("Accuracy ..." + ((double) (sumStatistics) / (double) countTimes));
    }

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    public static void CollapsedGibbsSamplingforDGMM(double[] Samples, int[] Class) {
        System.out.println("CollapsedGibbsSamplingforDGMM");
        double sumStatistics = 0;
        int countTimes = 0;
        ArrayList<ArrayList<Integer>> cacheZs = new ArrayList<>(); //an arraylist of arraylist with sample id and the list of all assined cluster
        for(int i = 0; i < Samples.length; i++)
        {
            cacheZs.add(new ArrayList<>());
        }
        
        Random r = new Random();
        ArrayList<ClusterDataStructure> list = new ArrayList<>();
        int maximumcluster = 1000;//just a stupid number
        for(int i = 0; i < maximumcluster; i++) {
            list.add(null);
        }
        list.set(0, new ClusterDataStructure(Samples.length));
        list.get(0).index = 0;        
        Random mr = new Random();
        double sum = 0;
        for(int i = 0; i < Samples.length; i++) {
            list.get(0).listID[i] = true;
            sum+=Samples[i];
        }
        list.get(0).count = Samples.length;
        //list.get(0).theta = r.nextDouble(); //in principle we don't need to specify this?!
        list.get(0).theta = sum/Samples.length; //in principle we don't need to specify this?!
        
        
        int[] sequencesZ = new int[Samples.length];
        //at the begining, all sequences Z sample to cluster 0
        
        
        double maximumiteration = 3000;
        for (int iteration = 0; iteration < maximumiteration; iteration++) {
            //sampling Z
            sequencesZ = CollapsedGibbsdiscreteSamplesDPMM(Samples, list, sequencesZ, maximumcluster);
            if (iteration > burnin) {
                if ((iteration - burnin) % 1==0) {                    
                    for (int i = 0; i < sequencesZ.length; i++) {
                        cacheZs.get(i).add(sequencesZ[i]);
                    }
                }
            
            if(iteration==(maximumiteration-1) || (iteration - burnin) % 100==0) {                
                double truecase = 0;
                int[] myclass = new int[Samples.length];
                for (int i = 0; i < sequencesZ.length; i++) {
                    
                    int[] counts = new int[maximumcluster];
                    if (cacheZs.get(i).size() > 100) { //not a great solution, I know
                        for (int j = cacheZs.get(i).size()-100; j < cacheZs.get(i).size(); j++) {
                            counts[cacheZs.get(i).get(j)] = counts[cacheZs.get(i).get(j)] + 1;
                        }
                    } else {
                        for (int j = 0; j < cacheZs.get(i).size(); j++) {
                            counts[cacheZs.get(i).get(j)] = counts[cacheZs.get(i).get(j)] + 1;
                        }
                    }
                    
                    
                    double dmax = -1.0;
                    int imax = -1;
                    for (int j = 0; j < maximumcluster; j++) {
                        double tempvalue = counts[j];
                        if (dmax < tempvalue) {
                            dmax = tempvalue;
                            imax = j;
                        }
                    }
                    myclass[i] = imax;
                }
                truecase = measure(Class, myclass);
                countTimes++;
                sumStatistics+=truecase;
                //System.out.println("Iteration " + iteration + "\nAccuracy ..." + ((double) (truecase)));
                int activeclass = 0;
                for (int i = 0; i < maximumcluster; i++) {
                    if (list.get(i) == null) {
                        continue;
                    }
                    activeclass++;

                }
                System.out.println("Iteration " + iteration + " Active classes" + activeclass);
                }
            }

            //System.out.println("Iteration "+iteration+" Active classes" + activeclass);
            
        }
        System.out.println("Accuracy ..." + ((double) (sumStatistics) / (double) countTimes));

    }

    
    
    
    
    
    
    public static void CollapsedGibbsSamplingforDGMMwithpartialSupervision(double[] Samples, int[] Class) {
        double sumStatistics = 0;
        int countTimes = 0;
        System.out.println("CollapsedGibbsSamplingforDGMMwithpartialSupervision");
        ArrayList<ArrayList<Integer>> cacheZs = new ArrayList<>(); //an arraylist of arraylist with sample id and the list of all assined cluster
        for(int i = 0; i < Samples.length; i++)
        {
            cacheZs.add(new ArrayList<>());
        }
        
        Random r = new Random();
        ArrayList<ClusterDataStructure> list = new ArrayList<>();
        int maximumcluster = 1000;//just a stupid number
        for(int i = 0; i < maximumcluster; i++) {
            list.add(null);
        }
        list.set(0, new ClusterDataStructure(Samples.length));
        list.get(0).index = 0;        
        Random mr = new Random();
        double sum = 0;
        for(int i = 0; i < Samples.length; i++) {
            list.get(0).listID[i] = true;
            sum+=Samples[i];
        }
        list.get(0).count = Samples.length;
        //list.get(0).theta = r.nextDouble(); //in principle we don't need to specify this?!
        list.get(0).theta = sum/Samples.length; //in principle we don't need to specify this?!
        
        
        int[] sequencesZ = new int[Samples.length];
        //at the begining, all sequences Z sample to cluster 0
        
        
        double maximumiteration = 3000;
        for (int iteration = 0; iteration < maximumiteration; iteration++) {
            //sampling Z
            sequencesZ = CollapsedGibbsdiscreteSamplesDPMMwithpartialSupervision(Samples, list, sequencesZ, maximumcluster, Class);
            if (iteration > burnin) {
                if ((iteration - burnin) % 1 == 0) {
                    for (int i = 0; i < sequencesZ.length; i++) {
                        cacheZs.get(i).add(sequencesZ[i]);
                    }
                }

                if (iteration == (maximumiteration - 1) || (iteration - burnin) % 100 == 0) {
                    double truecase = 0;
                    int[] myclass = new int[Samples.length];
                    for (int i = 0; i < sequencesZ.length; i++) {

                        int[] counts = new int[maximumcluster];
                        if (cacheZs.get(i).size() > 100) { //not a great solution, I know
                            for (int j = cacheZs.get(i).size() - 100; j < cacheZs.get(i).size(); j++) {
                                counts[cacheZs.get(i).get(j)] = counts[cacheZs.get(i).get(j)] + 1;
                            }
                        } else {
                            for (int j = 0; j < cacheZs.get(i).size(); j++) {
                                counts[cacheZs.get(i).get(j)] = counts[cacheZs.get(i).get(j)] + 1;
                            }
                        }
                        double dmax = -1.0;
                        int imax = -1;
                        for (int j = 0; j < maximumcluster; j++) {
                            double tempvalue = counts[j];
                            if (dmax < tempvalue) {
                                dmax = tempvalue;
                                imax = j;
                            }
                        }
                        myclass[i] = imax;
                    }
                    truecase = measure(Class, myclass);
                    sumStatistics+=truecase;
                    countTimes++;
                    
                    int activeclass = 0;
                    for (int i = 0; i < maximumcluster; i++) {
                        if (list.get(i) == null) {
                            continue;
                        }
                        activeclass++;

                    }
                    System.out.println("Iteration " + iteration + " Active classes" + activeclass);
                }
            }

            //System.out.println("Iteration "+iteration+" Active classes" + activeclass);
        }
        System.out.println("Accuracy ..." + ((double) (sumStatistics)/(double) countTimes));
        

    }

    
    public static int[] CollapsedGibbsdiscreteSamplesDPMM(double[] Samples, ArrayList<ClusterDataStructure> list, 
            int[] originalSequencesZ, int totalLimitClustersToDate) {
        int[] sequencesZ = new int[Samples.length];
        Random r = new Random();
        for (int sample = 0; sample < Samples.length; sample++) {
            double As[] = new double[totalLimitClustersToDate];
            int k = originalSequencesZ[sample];
            list.get(k).listID[sample] = false; //I remove it to do sample a new cluster.
            list.get(k).count = list.get(k).count - 1;
            if(list.get(k).count==0) {
                //ok I remove it
                list.set(k, null);
            }
            double[] Probs = new double[totalLimitClustersToDate];
            boolean nullforthefirstime = false;
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                if (list.get(i) != null) {

                    double mu_0 = 0.0;
                    double variance_0 = 1000;

                    double sum = 0;
                    for (int j = 0; j < list.get(i).listID.length; j++) {
                        if (list.get(i).listID[j] == true) {
                            sum += Samples[j];
                        }
                    }

                    double newMeansforsampling = varianceFIX * mu_0 + variance_0 * sum;
                    newMeansforsampling /= (varianceFIX + variance_0 * list.get(i).count);

                    //I ignore this because I assume mu_0 = 0 
                    double newVarianceforsampling = varianceFIX * variance_0;
                    newVarianceforsampling /= (varianceFIX + variance_0 * list.get(i).count);
                    

                    Probs[i] = integralFormula(Samples[sample], varianceFIX, newMeansforsampling, newVarianceforsampling);
                    //Probs[i] = probabilitydensityfunct(list.get(i).theta, varianceFIX, Samples[sample]);
                    Probs[i] *= (list.get(i).count) / (Samples.length + 1.0 - 1.0);
                } else {
                    //if null for the first time; ok new cluster
                    if (nullforthefirstime == false) {
                        nullforthefirstime = true; //i do this onlly once.
                        Probs[i] = integralFormula(Samples[sample], varianceFIX, 0, 100);
                        Probs[i] *= 1.0 / (Samples.length + 1.0 - 1.0);
                    }
                }
            }
            //normalizeProbs first
            double sum = 0;
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                sum += Probs[i];
            }
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                Probs[i] /= sum;
            }
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                if (i > 0) {
                    As[i] = As[i - 1] + Probs[i];
                } else {
                    As[i] = Probs[i];
                }
            }
            double rndNumber = r.nextDouble();
            LOOP:
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                if (i == 0) {
                    if (rndNumber <= As[i]) {
                        //pick i
                        sequencesZ[sample] = i;
                        if(list.get(i)==null) {
                            list.set(i, new ClusterDataStructure(Samples.length));
                            list.get(i).index = i;
                            //sample new mean - refer to jordan - need to proof

                            double newmean = varianceFIX * 0 + 100 * Samples[sample];
                            newmean /= (varianceFIX + 100);

                            //I ignore this because I assume mu_0 = 0 
                            double newvariance = varianceFIX * 100;
                            newvariance /= (varianceFIX + 100);
                            list.get(i).theta = getGaussian(newmean, Math.sqrt(newvariance));
                        }
                        list.get(i).listID[sample] = true;
                        list.get(i).count = list.get(i).count +1;
                        break LOOP;
                    }
                } else {
                    if (rndNumber > As[i - 1] && rndNumber <= As[i]) {            
                        if(i>25) {
                            System.out.println("Probably shit happens ...");
                        }
                        sequencesZ[sample] = i;
                        if(list.get(i)==null) {
                            list.set(i, new ClusterDataStructure(Samples.length));
                            list.get(i).index = i;
                            
                            double newmean = varianceFIX * 0 + 100 * Samples[sample];
                            newmean /= (varianceFIX + 100);

                            //sample new mean - refer to jordan - need to proof
                            double newvariance = varianceFIX * 100;
                            newvariance /= (varianceFIX + 100);
                            list.get(i).theta = getGaussian(newmean, Math.sqrt(newvariance));
                        }
                        list.get(i).listID[sample] = true;
                        list.get(i).count = list.get(i).count +1;
                        break LOOP;
                    }
                }
            }
        }
        return sequencesZ;
    }
    
    
    
    public static int[] CollapsedGibbsdiscreteSamplesDPMMwithpartialSupervision(double[] Samples, ArrayList<ClusterDataStructure> list, 
            int[] originalSequencesZ, int totalLimitClustersToDate, int[] Class) {
        int[] sequencesZ = new int[Samples.length];
        Random r = new Random();
        for (int sample = 0; sample < Samples.length; sample++) {
            double As[] = new double[totalLimitClustersToDate];
            int k = originalSequencesZ[sample];
            list.get(k).listID[sample] = false; //I remove it to do sample a new cluster.
            list.get(k).count = list.get(k).count - 1;
            if(list.get(k).count==0) {
                //ok I remove it
                list.set(k, null);
            }
            double[] Probs = new double[totalLimitClustersToDate];
            boolean nullforthefirstime = false;
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                if (list.get(i) != null) {

                    double mu_0 = 0.0;
                    double variance_0 = 1000;

                    double sum = 0;
                    for (int j = 0; j < list.get(i).listID.length; j++) {
                        if (list.get(i).listID[j] == true) {
                            sum += Samples[j];
                        }
                    }

                    double newMeansforsampling = varianceFIX * mu_0 + variance_0 * sum;
                    newMeansforsampling /= (varianceFIX + variance_0 * list.get(i).count);

                    //I ignore this because I assume mu_0 = 0 
                    double newVarianceforsampling = varianceFIX * variance_0;
                    newVarianceforsampling /= (varianceFIX + variance_0 * list.get(i).count);
                    

                    Probs[i] = integralFormula(Samples[sample], varianceFIX, newMeansforsampling, newVarianceforsampling);
                    //Probs[i] = probabilitydensityfunct(list.get(i).theta, varianceFIX, Samples[sample]);
                    Probs[i] *= (list.get(i).count) / (Samples.length + 1.0 - 1.0);
                } else {
                    //if null for the first time; ok new cluster
                    if (nullforthefirstime == false) {
                        nullforthefirstime = true; //i do this onlly once.
                        Probs[i] = integralFormula(Samples[sample], varianceFIX, 0, 100);
                        Probs[i] *= 1.0 / (Samples.length + 1.0 - 1.0);
                    }
                }
            }
            //normalizeProbs first
            double sum = 0;
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                sum += Probs[i];
            }
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                Probs[i] /= sum;
            }
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                if (i > 0) {
                    As[i] = As[i - 1] + Probs[i];
                } else {
                    As[i] = Probs[i];
                }
            }
            double rndNumber = r.nextDouble();
            LOOP:
            for (int i = 0; i < totalLimitClustersToDate; i++) {
                //supervised here
                if(Class[sample]==0) {
                    //class 1
                    if(sample%3==0) {
                        //one third belongs to class 0
                        //assign 0
                        i = Class[sample];
                        sequencesZ[sample] = i;
                        if(list.get(i)==null) {
                            list.set(i, new ClusterDataStructure(Samples.length));
                            list.get(i).index = i;
                            //sample new mean - refer to jordan - need to proof

                            double newmean = varianceFIX * 0 + 100 * Samples[sample];
                            newmean /= (varianceFIX + 100);

                            //I ignore this because I assume mu_0 = 0 
                            double newvariance = varianceFIX * 100;
                            newvariance /= (varianceFIX + 100);
                            list.get(i).theta = getGaussian(newmean, Math.sqrt(newvariance));
                        }
                        list.get(i).listID[sample] = true;
                        list.get(i).count = list.get(i).count +1;
                        break LOOP;
                        
                    }
                }
                
                
                if(Class[sample]==1) {
                    //class 1
                    if(sample%3==0) {
                        //one third belongs to class 0
                        //assign 0
                        i = Class[sample];
                        sequencesZ[sample] = i;
                        if(list.get(i)==null) {
                            list.set(i, new ClusterDataStructure(Samples.length));
                            list.get(i).index = i;
                            //sample new mean - refer to jordan - need to proof

                            double newmean = varianceFIX * 0 + 100 * Samples[sample];
                            newmean /= (varianceFIX + 100);

                            //I ignore this because I assume mu_0 = 0 
                            double newvariance = varianceFIX * 100;
                            newvariance /= (varianceFIX + 100);
                            list.get(i).theta = getGaussian(newmean, Math.sqrt(newvariance));
                        }
                        list.get(i).listID[sample] = true;
                        list.get(i).count = list.get(i).count +1;
                        break LOOP;
                        
                    }
                }
                if (i == 0) {
                    if (rndNumber <= As[i]) {
                        //pick i
                        sequencesZ[sample] = i;
                        if(list.get(i)==null) {
                            list.set(i, new ClusterDataStructure(Samples.length));
                            list.get(i).index = i;
                            //sample new mean - refer to jordan - need to proof

                            double newmean = varianceFIX * 0 + 100 * Samples[sample];
                            newmean /= (varianceFIX + 100);

                            //I ignore this because I assume mu_0 = 0 
                            double newvariance = varianceFIX * 100;
                            newvariance /= (varianceFIX + 100);
                            list.get(i).theta = getGaussian(newmean, Math.sqrt(newvariance));
                        }
                        list.get(i).listID[sample] = true;
                        list.get(i).count = list.get(i).count +1;
                        break LOOP;
                    }
                } else {
                    if (rndNumber > As[i - 1] && rndNumber <= As[i]) {            
                        sequencesZ[sample] = i;
                        if(list.get(i)==null) {
                            list.set(i, new ClusterDataStructure(Samples.length));
                            list.get(i).index = i;
                            
                            double newmean = varianceFIX * 0 + 100 * Samples[sample];
                            newmean /= (varianceFIX + 100);

                            //sample new mean - refer to jordan - need to proof
                            double newvariance = varianceFIX * 100;
                            newvariance /= (varianceFIX + 100);
                            list.get(i).theta = getGaussian(newmean, Math.sqrt(newvariance));
                        }
                        list.get(i).listID[sample] = true;
                        list.get(i).count = list.get(i).count +1;
                        break LOOP;
                    }
                }
            }
        }
        return sequencesZ;
    }
    
    public static void writeFile(double[] Samples, int[] Class) throws IOException {
        FileWriter samples = new FileWriter("Samples");
        for(int i = 0; i < Samples.length; i++) {
            samples.write(Samples[i]+"\t"+Class[i]+"\n");
        }       
        samples.close();
    }
    public static void loadFile(int totalPoints) throws FileNotFoundException, IOException {
        BufferedReader buf = new BufferedReader(new FileReader("Samples"));
        String s = "";
        Samples = new double[totalPoints];
        Class = new int[totalPoints];
        int count = 0;
        while((s=buf.readLine())!=null) {
            int k = s.indexOf("\t");
            Samples[count] = Double.parseDouble(s.substring(0, k).trim());
            Class[count] = Integer.parseInt(s.substring(k+1).trim());
            count++;
        }
        buf.close();
    }
    
    
    public static int[] unique(int[] x) {
        HashSet<Integer> hash = new HashSet<>();
        for (int i = 0; i < x.length; i++) {
            hash.add(x[i]);
        }

        int[] y = new int[hash.size()];

        Iterator<Integer> keys = hash.iterator();
        for (int i = 0; i < y.length; i++) {
            y[i] = keys.next();
        }

        return y;
    }
    
    
        /**
     * n choose k. Returns 0 if n is less than k.
     */
    public static double choose(int n, int k) {
        if (n < 0 || k < 0) {
            throw new IllegalArgumentException(String.format("Invalid n = %d, k = %d", n, k));
        }

        if (n < k) {
            return 0.0;
        }

        return Math.floor(0.5 + Math.exp(logChoose(n, k)));
    }
    
    /**
     * log of n choose k
     */
    public static double logChoose(int n, int k) {
        if (n < 0 || k < 0 || k > n) {
            throw new IllegalArgumentException(String.format("Invalid n = %d, k = %d", n, k));
        }

        return logFactorial(n) - logFactorial(k) - logFactorial(n - k);
    }
    
     /**
     * log of factorial of n
     */
    public static double logFactorial(int n) {
        if (n < 0) {
            throw new IllegalArgumentException(String.format("n has to be nonnegative: %d", n));
        }

        double f = 0.0;
        for (int i = 2; i <= n; i++) {
            f += Math.log(i);
        }

        return f;
    }
    
    public static double measure(int[] y1, int[] y2) {
        if (y1.length != y2.length) {
            throw new IllegalArgumentException(String.format("The vector sizes don't match: %d != %d.", y1.length, y2.length));
        }

        // Get # of non-zero classes in each solution
        int n = y1.length;

        int[] label1 = unique(y1);
        int n1 = label1.length;

        int[] label2 = unique(y2);
        int n2 = label2.length;

        // Calculate N contingency matrix
        int[][] count = new int[n1][n2];
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                int match = 0;

                for (int k = 0; k < n; k++) {
                    if (y1[k] == label1[i] && y2[k] == label2[j]) {
                        match++;
                    }
                }

                count[i][j] = match;
            }
        }

        // Marginals
        int[] count1 = new int[n1];
        int[] count2 = new int[n2];

        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                count1[i] += count[i][j];
                count2[j] += count[i][j];
            }
        }

        // Calculate RAND - Adj
        double rand1 = 0.0;
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                if (count[i][j] >= 2) {
                    rand1 += choose(count[i][j], 2);
                }
            }
        }

        double rand2a = 0.0;
        for (int i = 0; i < n1; i++) {
            if (count1[i] >= 2) {
                rand2a += choose(count1[i], 2);
            }
        }

        double rand2b = 0;
        for (int j = 0; j < n2; j++) {
            if (count2[j] >= 2) {
                rand2b += choose(count2[j], 2);
            }
        }

        double rand3 = rand2a * rand2b;
        rand3 /= choose(n, 2);
        double rand_N = rand1 - rand3;

        // D
        double rand4 = (rand2a + rand2b) / 2;
        double rand_D = rand4 - rand3;

        double rand = rand_N / rand_D;
        return rand;
    }

    public static double integralFormula(double x_i, double variance_fix, double mu_0, double variance_0) {
        double a = 0.5 / variance_0 + 0.5 / variance_fix;
        double b = mu_0 / variance_0 + x_i / variance_fix;
        double c = -1.0 * (mu_0 * mu_0 * 0.5 / variance_0 + x_i * x_i * 0.5 / variance_fix);
        double result = ((0.5 / Math.PI) / Math.sqrt(variance_0 * variance_fix)) * Math.sqrt(Math.PI / a) * Math.exp(b * b / 4.0 / a + c);
        //System.out.println(result);
        return result;
    }

    
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        //setting G_0 = N(0, 1)
        
        
        int[] a = new int[300];
        int[] b = new int[300];
        /*Random r = new Random();
        for (int i = 0; i < a.length; i++) {
            if(i%number_of_clusters==0) {
                a[i] = r.nextInt(number_of_clusters);
                b[i] = r.nextInt(number_of_clusters);
            }
            else {
                a[i] = 0;
                b[i] = 4;
            }

        }*/
        
        //System.out.println(measure(a, b));
        Samples = generateSamplesforToyProblems(8);
        writeFile(Samples, Class);
        loadFile(8 * 50);
        CollapsedGibbsSamplingforDGMMwithpartialSupervision(Samples, Class);
        CollapsedGibbsSamplingforDGMM(Samples, Class);
        GibbsSamplingforDGMM(Samples, Class);
        //GibbsSamplingforFMM(number_of_clusters, Samples, Class);
        //CollapsedGibbsSamplingforFMM(number_of_clusters, Samples, Class);
        //EM(number_of_clusters, Samples, Class);

    }

}


class Statistics 
{
    double[] data;
    int size;   

    public Statistics(double[] data) 
    {
        this.data = data;
        size = data.length;
    }   

    double getMean()
    {
        double sum = 0.0;
        for(double a : data)
            sum += a;
        return sum/size;
    }

    double getVariance()
    {
        double mean = getMean();
        double temp = 0;
        for(double a :data)
            temp += (a-mean)*(a-mean);
        return temp/size;
    }

    double getStdDev()
    {
        return Math.sqrt(getVariance());
    }

    public double median() 
    {
       Arrays.sort(data);

       if (data.length % 2 == 0) 
       {
          return (data[(data.length / 2) - 1] + data[data.length / 2]) / 2.0;
       } 
       return data[data.length / 2];
    }
}