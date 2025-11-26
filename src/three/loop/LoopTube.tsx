import { Tube } from '@react-three/drei'
import { useMemo } from 'react'
import { CatmullRomCurve3, Vector3 } from 'three'
import { EARTH_RENDER_RADIUS, LOS_ANGELES, SAN_FRANCISCO, interpolateGreatCircle } from '../../config/simulation'

const LoopTube = () => {
  const curve = useMemo(() => {
    const points: Vector3[] = []
    const segments = 24
    for (let i = 0; i <= segments; i += 1) {
      points.push(interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, i / segments, EARTH_RENDER_RADIUS + 0.02))
    }
    return new CatmullRomCurve3(points)
  }, [])

  return (
    <Tube args={[curve, 128, 0.01, 12, false]}>
      <meshStandardMaterial color="#6aa0ff" transparent opacity={0.4} />
    </Tube>
  )
}

export default LoopTube
