using UnityEngine;

namespace Avionic.Math
{
    public class UniverseContainer : MonoBehaviour
    {
        public double EditorScale = 1/10000D;
        
        public Vector3d Offset
        {
            get => _offset;
            set
            {
                _offset = value;
                OffsetChanged?.Invoke(this);
            }
        }

        public double Scale
        {
            get => _scale;
            set
            {
                _scale = value;
                ScaleChanged?.Invoke(this);
            }
        }

        private double _scale;
        private Vector3d _offset;
        
        public delegate void Callback(UniverseContainer container);
        public event Callback ScaleChanged;
        public event Callback OffsetChanged;

        private void Awake()
        {
            Scale = EditorScale;
        }
    }
}