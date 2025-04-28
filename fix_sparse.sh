for f in *; do
  if [ -f "$f" ]; then
    cp --sparse=never "$f" "$f.tmp" && mv "$f.tmp" "$f"
  fi
done
