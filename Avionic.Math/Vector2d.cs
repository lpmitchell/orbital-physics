using System;
using System.Diagnostics.CodeAnalysis;
using UnityEngine;

namespace Avionic.Math
{
    [SuppressMessage("ReSharper", "InconsistentNaming"), Serializable]
    public struct Vector2d
    {
        public double x;
        public double y;
        
        private const double KEpsilon = 0.0000000001;
        
        public double magnitude => System.Math.Sqrt(x * x + y * y);
        public double sqrMagnitude => x * x + y * y;

        public Vector2d(double x, double y)
        {
            this.x = x;
            this.y = y;
        }
        
        public Vector2d normalized {
            get
            {
                var mag = magnitude;
                if (mag > KEpsilon)
                    return this / mag;
                return new Vector2d(0,0);
            }
        }
    
        public static double Dot(Vector2d vector1, Vector2d vector2)
        {
            return vector1.x * vector2.x +
                   vector1.y * vector2.y;
        }
        
        public static Vector2d operator+(Vector2d a, Vector2d b) => new(a.x + b.x, a.y + b.y);
        public static Vector2d operator-(Vector2d a, Vector2d b) => new(a.x - b.x, a.y - b.y);
        public static Vector2d operator-(Vector2d a)             => new(-a.x, -a.y);
        public static Vector2d operator*(Vector2d a, double d)   => new(a.x * d, a.y * d);
        public static Vector2d operator*(double d, Vector2d a)   => new(a.x * d, a.y * d);
        public static Vector2d operator*(Vector2d a, Vector2d b) => new(a.x * b.x, a.y * b.y);
        public static Vector2d operator/(Vector2d a, double d)   => new(a.x / d, a.y / d);

        public static Vector2d midpoint(Vector2d a, Vector2d b) => new((a.x + b.x) / 2d, (a.y + b.y) / 2d);
        public static Vector2d midpoint(Vector2d a, Vector2d b, Vector2d c) => new((a.x + b.x + c.x) / 3d, (a.y + b.y + c.y) / 3d);
        public static Vector2d midpoint(Vector2d a, Vector2d b, Vector2d c, Vector2d d) => new((a.x + b.x + c.x + d.x) / 4d, (a.y + b.y + c.y + d.y) / 4d);

        public void toMidpoint(Vector2d b)
        {
            x += b.x;
            y += b.y;
            x /= 2d;
            y /= 2d;
        }

        public static implicit operator Vector2(Vector2d t) => new((float)t.x, (float)t.y);
        public static implicit operator Vector2d(Vector2 t) => new(t.x, t.y);
        public static implicit operator Vector3(Vector2d t) => new((float)t.x, (float)t.y, 0f);

        public override string ToString() => $"Vector2d({x:F4}, {y:F4})";
    
        [SuppressMessage("ReSharper", "NonReadonlyMemberInGetHashCode")]
        public override int GetHashCode()
        {
            return (int)(x * 1000000) ^ (int)(y * 1000000) << 8;
        }

        public static Vector2d Lerp(Vector2d a, Vector2d b, double t)
        {
            if (t < 0) t = 0;
            if (t > 1) t = 1;
            return new Vector2d(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
        }
        
        public (double, double) IntersectTriangleUvEpsilon(Vector2 a, Vector2 b, Vector2 c, double negativeEpsilon = -MathConstants.Epsilon)
        {
            var detT = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
            var alpha = ((b.y - c.y) * (x - c.x) + (c.x - b.x) * (y - c.y)) / detT;
            var beta = ((c.y - a.y) * (x - c.x) + (a.x - c.x) * (y - c.y)) / detT;
            var gamma = 1 - alpha - beta;

            if (alpha >= negativeEpsilon && beta >= negativeEpsilon && gamma >= negativeEpsilon)
            {
                return (alpha, beta);
            }

            return (double.NaN, double.NaN);
        }
        
        public (double, double) IntersectTriangleUv(Vector2 a, Vector2 b, Vector2 c)
        {
            var detT = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
            var alpha = ((b.y - c.y) * (x - c.x) + (c.x - b.x) * (y - c.y)) / detT;
            var beta = ((c.y - a.y) * (x - c.x) + (a.x - c.x) * (y - c.y)) / detT;
            var gamma = 1 - alpha - beta;

            if (alpha >= 0 && beta >= 0 && gamma >= 0)
            {
                return (alpha, beta);
            }

            return (double.NaN, double.NaN);
        }
        
        public (Vector2d, double) ClosestPointInTriangle(Vector2d a, Vector2d b, Vector2d c)
        {
            var abClosest = ClosestPointOnLineSegment(a, b);
            var bcClosest = ClosestPointOnLineSegment(b, c);
            var caClosest = ClosestPointOnLineSegment(c, a);

            var dAB = (abClosest - this).sqrMagnitude;
            var dBC = (bcClosest - this).sqrMagnitude;
            var dCA = (caClosest - this).sqrMagnitude;

            if (dAB <= dBC && dAB <= dCA) return (abClosest, dAB);
            if (dBC <= dAB && dBC <= dCA) return (bcClosest, dBC);
            return (caClosest, dCA);
        }
        
        public Vector2d ClosestPointOnLineSegment(Vector2d p1, Vector2d p2)
        {
            var p1ToP2 = p2 - p1;
            var p1ToPoint = this - p1;

            var t = Dot(p1ToPoint, p1ToP2) / Dot(p1ToP2, p1ToP2);
            t = System.Math.Max(0, System.Math.Min(1, t));

            return p1 + t * p1ToP2;
        }

    }
}