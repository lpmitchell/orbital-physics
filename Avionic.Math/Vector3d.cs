using System;
using System.Diagnostics.CodeAnalysis;
using Unity.Mathematics;
using UnityEngine;

namespace Avionic.Math
{
    [SuppressMessage("ReSharper", "InconsistentNaming"), Serializable]
    public struct Vector3d
    {
        public double x;
        public double y;
        public double z;
        
        private const double KEpsilon = 0.0000000001;
        
        public double magnitude => System.Math.Sqrt(x * x + y * y + z * z);
        public double sqrMagnitude => x * x + y * y + z * z;

        public Vector3d(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
        
        public Vector3d normalized {
            get
            {
                var mag = magnitude;
                if (mag > KEpsilon)
                    return this / mag;
                return new Vector3d(0,0,0);
            }
        }
        
        public static double Angle(Vector3d from, Vector3d to)
        {
            var dot = Dot(from.normalized, to.normalized);
            return System.Math.Acos(dot < -1.0 ? -1.0 : (dot > 1.0 ? 1.0 : dot)) * 57.29578d;
        }
    
        public static Vector3d Cross(Vector3d vector1, Vector3d vector2)
        {
            return new Vector3d(
                vector1.y * vector2.z - vector1.z * vector2.y,
                vector1.z * vector2.x - vector1.x * vector2.z,
                vector1.x * vector2.y - vector1.y * vector2.x);
        }
    
        public static double Dot(Vector3d vector1, Vector3d vector2)
        {
            return vector1.x * vector2.x +
                   vector1.y * vector2.y +
                   vector1.z * vector2.z;
        }
        
        public static Vector3d operator+(Vector3d a, Vector3d b) => new(a.x + b.x, a.y + b.y, a.z + b.z);
        public static Vector3d operator-(Vector3d a, Vector3d b) => new(a.x - b.x, a.y - b.y, a.z - b.z);
        public static Vector3d operator-(Vector3d a)             => new(-a.x, -a.y, -a.z);
        public static Vector3d operator*(Vector3d a, double d)   => new(a.x * d, a.y * d, a.z * d);
        public static Vector3d operator*(double d, Vector3d a)   => new(a.x * d, a.y * d, a.z * d);
        public static Vector3d operator/(Vector3d a, double d)   => new(a.x / d, a.y / d, a.z / d);
        
        
        public static Vector3d operator+(Vector3d a, Vector3 b) => new(a.x + b.x, a.y + b.y, a.z + b.z);
        public static Vector3d operator-(Vector3d a, Vector3 b) => new(a.x - b.x, a.y - b.y, a.z - b.z);
        public static Vector3d operator+(Vector3 a, Vector3d b) => new(a.x + b.x, a.y + b.y, a.z + b.z);
        public static Vector3d operator-(Vector3 a, Vector3d b) => new(a.x - b.x, a.y - b.y, a.z - b.z);

        public static Vector3d midpoint(Vector3d a, Vector3d b) => new((a.x + b.x) / 2d, (a.y + b.y) / 2d, (a.z + b.z) / 2d);
        public static Vector3d midpoint(Vector3d a, Vector3d b, Vector3d c) => new((a.x + b.x + c.x) / 3d, (a.y + b.y + c.y) / 3d, (a.z + b.z + c.z) / 3d);
        public static Vector3d midpoint(Vector3d a, Vector3d b, Vector3d c, Vector3d d) => new((a.x + b.x + c.x + d.x) / 4d, (a.y + b.y + c.y + d.y) / 4d, (a.z + b.z + c.z + d.z) / 4d);

        public void toMidpoint(Vector3d b)
        {
            x += b.x;
            y += b.y;
            z += b.z;
            x /= 2d;
            y /= 2d;
            z /= 2d;
        }

        public bool IsWithinEpsilon(Vector3d other)
        {
            return (this - other).sqrMagnitude < MathConstants.Epsilon;
        }

        public static implicit operator Vector3(Vector3d t) => new((float)t.x, (float)t.y, (float)t.z);
        public static implicit operator Vector3d(Vector3 t) => new(t.x, t.y, t.z);
        public static implicit operator float3(Vector3d t) => new((float)t.x, (float)t.y, (float)t.z);
        public static implicit operator Vector3d(float3 t) => new(t.x, t.y, t.z);
        public static implicit operator double3(Vector3d t) => new(t.x, t.y, t.z);
        public static implicit operator Vector3d(double3 t) => new(t.x, t.y, t.z);
        
        public override string ToString() => $"Vector3d({x:F4}, {y:F4}, {z:F4})";
    
        [SuppressMessage("ReSharper", "NonReadonlyMemberInGetHashCode")]
        public override int GetHashCode()
        {
            return (int)(x * 1000000) ^ (int)(y * 1000000) << 8 ^ (int)(z * 1000000) << 16;
        }

        public static Vector3d Lerp(Vector3d a, Vector3d b, double t)
        {
            if (t < 0) t = 0;
            if (t > 1) t = 1;
            return new Vector3d(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t);
        }
    }
}