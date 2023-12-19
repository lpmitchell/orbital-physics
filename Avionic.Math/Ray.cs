using System.Text;
using UnityEngine;

namespace Avionic.Math
{
    public struct Ray
    {
        private const double KEpsilon = 0.0000000001;

        public Vector3d Direction;
        public Vector3d Origin;

        /// <summary>
        /// Ray-versus-triangle intersection test suitable for ray-tracing etc.
        /// Port of Möller–Trumbore algorithm c++ version from:
        /// https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm
        /// </summary>
        public double IntersectTriangle(Vector3d v0, Vector3d v1, Vector3d v2) {
            var e1 = v1 - v0;
            var e2 = v2 - v0;

            var h = Vector3d.Cross(Direction, e2);
            var a = Vector3d.Dot(e1, h);
            if (a > -KEpsilon && a < KEpsilon)
            {
                return double.NaN;
            }

            var f = 1.0 / a;

            var s = Origin - v0;
            var u = f * Vector3d.Dot(s, h);
            if (u < 0.0 || u > 1.0)
            {
                return double.NaN;
            }

            var q = Vector3d.Cross(s, e1);
            var v = f * Vector3d.Dot(Direction, q);
            if (v < 0.0 || u + v > 1.0)
            {
                return double.NaN;
            }

            var t = f * Vector3d.Dot(e2, q);
            return t > KEpsilon ? t : double.NaN;
        }

        private static StringBuilder _log = new StringBuilder();
        
        public bool TriangleIsApproximatelyInDirection(Vector3d v0, Vector3d v1, Vector3d v2, double buffer)
        {
            var mid = Vector3d.midpoint(v0, v1, v2).normalized;
            var size = (v0.normalized - mid).sqrMagnitude;
            var projected = Direction.normalized;
            var dist = (projected - mid).sqrMagnitude;
            
            return dist / size < (2f + buffer);
        }

        public void DumpLog()
        {
            Debug.Log(_log);
            _log.Clear();
        }
        
        public double IntersectTriangleWithBuffer(Vector3d v0, Vector3d v1, Vector3d v2, double buffer = 0.1)
        {
            var e1 = v1 - v0;
            var e2 = v2 - v0;

            var h = Vector3d.Cross(Direction, e2);
            var a = Vector3d.Dot(e1, h);
            if (a > -KEpsilon && a < KEpsilon)
            {
                return double.NaN;
            }

            var f = 1.0 / a;

            var s = Origin - v0;
            var u = f * Vector3d.Dot(s, h);
            if (u < -buffer || u > 1.0 + buffer)
            {
                return double.NaN;
            }

            var q = Vector3d.Cross(s, e1);
            var v = f * Vector3d.Dot(Direction, q);
            if (v < -buffer || u + v > 1.0 + buffer)
            {
                return double.NaN;
            }

            var t = f * Vector3d.Dot(e2, q);
            return t > KEpsilon ? t : double.NaN;
        }
        
        public (double, double) IntersectTriangleUv(Vector3d v0, Vector3d v1, Vector3d v2) {
            var e1 = v1 - v0;
            var e2 = v2 - v0;

            var h = Vector3d.Cross(Direction, e2);
            var a = Vector3d.Dot(e1, h);
            if (a > -KEpsilon && a < KEpsilon)
            {
                return (double.NaN, double.NaN);
            }

            var f = 1.0 / a;

            var s = Origin - v0;
            var u = f * Vector3d.Dot(s, h);
            if (u < 0.0 || u > 1.0)
            {
                return (double.NaN, double.NaN);
            }

            var q = Vector3d.Cross(s, e1);
            var v = f * Vector3d.Dot(Direction, q);
            if (v < 0.0 || u + v > 1.0)
            {
                return (double.NaN, double.NaN);
            }

            var t = f * Vector3d.Dot(e2, q);
            return t > KEpsilon ? (u,v) : (double.NaN, double.NaN);
        }
        
        public (double, double) IntersectTriangleUvFast(Vector3d v0, Vector3d v1, Vector3d v2) {
            var e1 = v1 - v0;
            var e2 = v2 - v0;

            var h = Vector3d.Cross(Direction, e2);
            var a = Vector3d.Dot(e1, h);
            if (a > -KEpsilon && a < KEpsilon)
            {
                return (double.NaN, double.NaN);
            }

            var f = 1.0 / a;

            var s = -v0;
            var u = f * Vector3d.Dot(s, h);
            if (u < 0.0 || u > 1.0)
            {
                return (double.NaN, double.NaN);
            }

            var q = Vector3d.Cross(s, e1);
            var v = f * Vector3d.Dot(Direction, q);
            if (v < 0.0 || u + v > 1.0)
            {
                return (double.NaN, double.NaN);
            }

            var t = f * Vector3d.Dot(e2, q);
            return t > KEpsilon ? (u,v) : (double.NaN, double.NaN);
        }
        
        public double IntersectTriPlane(Vector3d v0, Vector3d v1, Vector3d v2) {
            var e1 = v1 - v0;
            var e2 = v2 - v0;

            var h = Vector3d.Cross(Direction, e2);
            var a = Vector3d.Dot(e1, h);
            if (a > -KEpsilon && a < KEpsilon)
            {
                return double.NaN;
            }

            var f = 1.0 / a;

            var s = Origin - v0;
            var q = Vector3d.Cross(s, e1);
            
            var t = f * Vector3d.Dot(e2, q);
            return t > KEpsilon ? t : double.NaN;
        }

        public Vector3d GetPoint(double distance)
        {
            return Origin + Direction * distance;
        }
    }
}