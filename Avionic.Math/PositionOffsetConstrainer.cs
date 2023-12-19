using UnityEngine;

namespace Avionic.Math
{
    public class PositionOffsetConstrainer : MonoBehaviour
    {
        public double MaxDistanceFromOrigin = 2500;
        public UniverseContainer Universe;
        public PositionD Position;

        private void Awake()
        {
            Universe.Offset = -Position.Value;
        }

        private void LateUpdate()
        {
            if (transform.position.sqrMagnitude < MaxDistanceFromOrigin * MaxDistanceFromOrigin) return;
            Universe.Offset = -Position.Value;
        }
    }
}