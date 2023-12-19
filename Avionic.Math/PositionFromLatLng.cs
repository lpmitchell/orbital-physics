using UnityEngine;

namespace Avionic.Math
{
    public class PositionFromLatLng : MonoBehaviour
    {
        public double Radius;
        public PositionD Position;
        public Vector2 LatLng;
        public bool Rotate;
        private const double Deg2Rad = System.Math.PI / 180;

        private void Start()
        {
            var pos = CoordinateToPoint(LatLng);
            Position.Value = pos * Radius;
            if (Rotate) Position.transform.up = pos;
        }
        
        private static Vector3d CoordinateToPoint(Vector2 coordinate)
        {
            var y = System.Math.Sin(coordinate.x * Deg2Rad);
            var r = System.Math.Cos(coordinate.x * Deg2Rad); // radius of 2d circle cut through sphere at 'y'
            var x = System.Math.Sin(coordinate.y * Deg2Rad) * r;
            var z = -System.Math.Cos(coordinate.y * Deg2Rad) * r;

            return new Vector3d(x, y, z);
        }
    }
}