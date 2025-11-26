import { Line } from '@react-three/drei'
import { useMemo } from 'react'
import { interpolateGreatCircle, LOS_ANGELES, SAN_FRANCISCO } from '../../config/simulation'

const LoopTethers = () => {
  const tetherLines = useMemo(() => {
    return new Array(6).fill(0).map((_, idx) => {
      const t = idx / 5
      const anchor = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, t, 1.05)
      const ground = anchor.clone().setLength(1.02)
      return [anchor, ground]
    })
  }, [])

  return (
    <group>
      {tetherLines.map((pair, idx) => (
        <Line key={idx} points={pair.map((p) => p.toArray())} color="#b8c0ff" lineWidth={1.5} />
      ))}
    </group>
  )
}

export default LoopTethers
