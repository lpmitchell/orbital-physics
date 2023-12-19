using UnityEngine;

namespace Avionic.Math
{
    public static class VectorArrayHelpers
    {
        // Function to convert a Vector3d array to a Vector3 array
        public static Vector3[] ToVector3(this Vector3d[] array)
        {
            var r = new Vector3[array.Length];
            for (var i = 0; i < array.Length; i++)
            {
                r[i] = array[i];
            }
            return r;
        }
        
        // Function to convert a Vector3d array to a Vector3 array with an offset
        public static Vector3[] ToVector3(this Vector3d[] array, Vector3d offset)
        {
            var r = new Vector3[array.Length];
            for (var i = 0; i < array.Length; i++)
            {
                r[i] = array[i] + offset;
            }
            return r;
        }
        
        // Function to convert a Vector3d array to a Vector3 array with an offset and scale
        public static Vector3[] ToVector3(this Vector3d[] array, Vector3d offset, double scale)
        {
            var r = new Vector3[array.Length];
            for (var i = 0; i < array.Length; i++)
            {
                r[i] = (array[i] + offset) * scale;
            }
            return r;
        }
        
        // Function to convert a Vector2d array to a Vector2 array
        public static Vector2[] ToVector2(this Vector2d[] array)
        {
            var r = new Vector2[array.Length];
            for (var i = 0; i < array.Length; i++)
            {
                r[i] = array[i];
            }
            return r;
        }
    }
}