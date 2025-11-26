import { Instances, Instance } from '@react-three/drei'
import { useMemo } from 'react'
import { interpolateGreatCircle, LOS_ANGELES, SAN_FRANCISCO, LOOP_RADIUS_KM } from '../../config/simulation'

const LoopSections = () => {
  const points = useMemo(() => {
    return new Array(12).fill(0).map((_, idx) => {
      const t = idx / 11
      return interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, t, LOOP_RADIUS_KM)
    })
  }, [])

  // Section spheres: 0.015 km = 15 meters radius (representing buoyant vacuum segments)
  return (
    <Instances limit={points.length} position={[0, 0, 0]}>
      <sphereGeometry args={[0.015, 12, 12]} />
      <meshStandardMaterial color="#7af2b5" emissive="#1b462f" emissiveIntensity={0.6} />
      {points.map((point, idx) => (
        <Instance key={idx} position={point.toArray()} />
      ))}
    </Instances>
  )
}

export default LoopSections
