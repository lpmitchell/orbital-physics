using Avionic.Math;
using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;

namespace Avionic.Simulation.Kepler
{
    [BurstCompile]
    public struct PositionSolverJobElliptical : IJobParallelFor
    {
        [ReadOnly] public double StartTime;
        [ReadOnly] public double MeanAnomaly;
        [ReadOnly] public double MeanMotion;
        [ReadOnly] public double Step;
        [ReadOnly] public double Eccentricity;
        [ReadOnly] public double SemiMinorAxis;
        [ReadOnly] public double SemiMajorAxis;
        [ReadOnly] public Vector3d SemiMinorAxisBasis;
        [ReadOnly] public Vector3d SemiMajorAxisBasis;
        [ReadOnly] public Vector3d Offset;
        
        [WriteOnly] public NativeArray<Vector3d> Result;

        private double SolveKepler(double m, double e)
        {
            var origAnomaly = m;
            var iterations = (int)math.ceil((e + 0.7) * 1.25) << 1;
            for (var i = 0; i < iterations; i++)
            {
                var esinE = e * math.sin(m);
                var ecosE = e * math.cos(m);
                var deltaE = m - esinE - origAnomaly;
                var n = 1.0 - ecosE;
                m += -5d * deltaE / (n + math.sign(n) * math.sqrt(math.abs(16d * n * n - 20d * deltaE * esinE)));
            }
            return m;
        }

        public void Execute(int index)
        {
            var m = MeanAnomaly + MeanMotion * (StartTime + Step * index);
            var eccentricAnomaly = SolveKepler(m, Eccentricity);
            var rX = math.sin(eccentricAnomaly) * SemiMinorAxis;
            var rY = -math.cos(eccentricAnomaly) * SemiMajorAxis;
            Result[index] = (-SemiMinorAxisBasis * rX - SemiMajorAxisBasis * rY) + Offset;
        }
    }
    
    [BurstCompile]
    public struct PositionSolverJobHyperbolic : IJobParallelFor
    {
        [ReadOnly] public double StartTime;
        [ReadOnly] public double MeanAnomaly;
        [ReadOnly] public double MeanMotion;
        [ReadOnly] public double Step;
        [ReadOnly] public double Eccentricity;
        [ReadOnly] public double SemiMinorAxis;
        [ReadOnly] public double SemiMajorAxis;
        [ReadOnly] public Vector3d SemiMinorAxisBasis;
        [ReadOnly] public Vector3d SemiMajorAxisBasis;
        [ReadOnly] public Vector3d Offset;
        
        [WriteOnly] public NativeArray<Vector3d> Result;

        private double SolveKeplerHyperbolic(double m, double e)
        {
            var delta = 1d;

            // Adjusting Danby's guess for negative mean anomalies.
            var signM = m >= 0 ? 1d : -1d;
            var f = math.log(2d * math.abs(m) / e + 1.8d) * signM;
            if (double.IsNaN(f) || double.IsInfinity(f))
            {
                return m;
            }

            while (math.abs(delta) > 1e-8)
            {
                delta = (e * math.sinh(f) - f - m) / (e * math.cosh(f) - 1d);
                f -= delta;
            }
            
            return f;
        }

        public void Execute(int index)
        {
            var m = MeanAnomaly + MeanMotion * (StartTime + Step * index);
            var eccentricAnomaly = SolveKeplerHyperbolic(m, Eccentricity);

            var rX = math.sinh(eccentricAnomaly) * SemiMinorAxis;
            var rY = math.cosh(eccentricAnomaly) * SemiMajorAxis;
            Result[index] = (-SemiMinorAxisBasis * rX - SemiMajorAxisBasis * rY) + Offset;
        }
    }
    
    
    
    [BurstCompile]
    public struct PositionSolverJobCircular : IJobParallelFor
    {
        [ReadOnly] public double StartTime;
        [ReadOnly] public double MeanAnomaly;
        [ReadOnly] public double MeanMotion;
        [ReadOnly] public double Step;
        [ReadOnly] public double Eccentricity;
        [ReadOnly] public double PeriapsisDistance;
        [ReadOnly] public Vector3d SemiMinorAxisBasis;
        [ReadOnly] public Vector3d SemiMajorAxisBasis;
        [ReadOnly] public Vector3d Offset;
        
        [WriteOnly] public NativeArray<Vector3d> Result;

        private double SolveKeplerCircular(double m, double e)
        {
            var origAnomaly = m;
            var iterations = (int)math.ceil((e + 0.7) * 1.25) << 1;
            for (var i = 0; i < iterations; i++)
            {
                var esinE = e * math.sin(m);
                var ecosE = e * math.cos(m);
                var deltaE = m - esinE - origAnomaly;
                var n = 1.0 - ecosE;
                m += -5d * deltaE / (n + math.sign(n) * math.sqrt(math.abs(16d * n * n - 20d * deltaE * esinE)));
            }
            return m;
        }

        public void Execute(int index)
        {
            var m = MeanAnomaly + MeanMotion * (StartTime + Step * index);
            var eccentricAnomaly = SolveKeplerCircular(m, Eccentricity);
            var rX = math.sin(eccentricAnomaly) * PeriapsisDistance;
            var rY = -math.cos(eccentricAnomaly) * PeriapsisDistance;
            Result[index] = (-SemiMinorAxisBasis * rX - SemiMajorAxisBasis * rY) + Offset;
        }
    }
}