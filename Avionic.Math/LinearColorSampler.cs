using System;
using Unity.Collections;
using UnityEngine;

namespace Avionic.Math
{
    public struct LinearColorSampler: IDisposable
    {
        [ReadOnly] private NativeArray<Color> _values;
        [ReadOnly] private NativeArray<float> _positions;
        
        public LinearColorSampler(NativeArray<Color> values, NativeArray<float> positions)
        {
            _values = values;
            _positions = positions;
        }

        public void Dispose()
        {
            _values.Dispose();
            _positions.Dispose();
        }
        
        public Color Sample(float position)
        {
            if (position <= _positions[0]) return _values[0];
            if (position >= _positions[_positions.Length - 1]) return _values[_values.Length - 1];
            
            var index = 1;
            for (var i = 1; i < _positions.Length; i++)
            {
                if (position >= _positions[i]) continue;
                index = i;
                break;
            }

            var t = (position - _positions[index - 1]) / (_positions[index] - _positions[index - 1]);
            return Color.Lerp(_values[index - 1], _values[index], t);
        }
    }
}