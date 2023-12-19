using UnityEngine;

namespace Avionic.Math
{
    public class PositionFromDirection : MonoBehaviour
    {
        public double Radius;
        public PositionD Position;
        public Vector3 Direction;
        public bool Rotate;

        private void Start()
        {
            Position.Value = (Vector3d)Direction.normalized * Radius;
            if (Rotate) Position.transform.up = Direction;
        }
    }
}